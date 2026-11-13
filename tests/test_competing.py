import copy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest
from encounters.schema import Dataset
from encounters.values import digest
from encounters.journal import Journal
from pathways.report import analyze, verify
from pathways.query import dispatch
from pathways.provenance import competing_lineage
from waiting_analysis.competing import (OUTCOMES, as_fraction, cumulative_incidence, target_at,
                                        verify as verify_curve, verify_study, aligned_naive,
                                        naive_estimates)
from waiting_analysis.plan import run as run_plan, verify as verify_plan


def micro(visits, people=('q1',), facilities=None, referrals=None):
    people_rows = [{'id': key, 'county': 'east', 'age_band': 'older', 'rural': False} for key in people]
    facilities = facilities or [{'id': 'a', 'county': 'east', 'stage': 'acute', 'name': 'acute site'}]
    return {'version': 1, 'facilities': facilities, 'people': people_rows,
            'visits': visits, 'referrals': referrals or []}


def visit(key, person, admitted, discharged=None, outcome=None, recorded=None, facility='a'):
    return {'id': key, 'person': person, 'facility': facility, 'admitted': admitted,
            'discharged': discharged, 'outcome': outcome,
            'recorded': recorded or (discharged or admitted)}


def episode(key, duration, outcome):
    return {'id': key, 'duration': duration, 'outcome': outcome}


class CompetingRiskEstimationTests(unittest.TestCase):
    def test_same_day_multiple_outcomes_conserve_probability(self):
        curve = cumulative_incidence([
            episode('a', 3, 'completed'), episode('b', 3, 'death'),
            episode('c', 3, 'transfer'), episode('d', 3, None),
            episode('e', 6, 'completed')])
        first = curve['rows'][0]
        self.assertEqual(first['at_risk'], 5)
        self.assertEqual(first['outcome_events'], {'completed': 1, 'transfer': 1, 'death': 1, 'left': 0})
        self.assertEqual(first['events_total'], 3)
        self.assertEqual(first['censored'], 1)
        self.assertEqual(as_fraction(first['survival']), Fraction(2, 5))
        self.assertEqual(as_fraction(first['cumulative']['completed']), Fraction(1, 5))
        self.assertEqual(as_fraction(first['cumulative']['death']), Fraction(1, 5))
        self.assertEqual(as_fraction(first['cumulative']['transfer']), Fraction(1, 5))
        self.assertEqual(as_fraction(first['conservation']), 1)
        self.assertNotEqual(first['survival'], first['survival_before'])
        for row in curve['rows']:
            total = as_fraction(row['survival']) + sum((as_fraction(value) for value in row['cumulative'].values()), Fraction(0))
            self.assertEqual(total, 1)
        self.assertTrue(verify_curve(curve)['valid'])
        self.assertEqual(curve['events_total'], 4)
        self.assertEqual(curve['censored'], 1)
        self.assertEqual(curve['last_time'], 6)

    def test_all_open_episodes_censor_without_inventing_completion(self):
        curve = cumulative_incidence([episode('a', 4, None), episode('b', 9, None)])
        self.assertEqual(len(curve['rows']), 2)
        for row in curve['rows']:
            self.assertEqual(as_fraction(row['survival']), 1)
            self.assertEqual(as_fraction(row['cumulative']['completed']), 0)
            self.assertTrue(row['event_time'] is False)
        self.assertTrue(target_at(curve, 9)['estimable'])
        target = target_at(curve, 10)
        self.assertFalse(target['estimable'])
        self.assertEqual(target['reason'], 'followup_ends_before_target')
        self.assertIsNone(target['cumulative'])
        self.assertEqual(target_at(curve, 2)['cumulative']['completed'], {'numerator': 0, 'denominator': 1})
        self.assertTrue(verify_curve(curve)['valid'])

    def test_only_competing_outcomes_leave_no_recovery_survival(self):
        curve = cumulative_incidence([episode('a', 2, 'death'), episode('b', 2, 'transfer'),
                                      episode('c', 5, 'left')])
        first = curve['rows'][0]
        self.assertEqual(as_fraction(first['survival']), Fraction(1, 3))
        self.assertEqual(as_fraction(first['cumulative']['completed']), 0)
        target = target_at(curve, 2)
        self.assertEqual(as_fraction(target['cumulative']['death']), Fraction(1, 3))
        self.assertEqual(as_fraction(target['cumulative']['transfer']), Fraction(1, 3))
        self.assertEqual(target['cumulative']['completed'], {'numerator': 0, 'denominator': 1})
        self.assertEqual(as_fraction(target['survival']), Fraction(1, 3))
        self.assertTrue(verify_curve(curve)['valid'])
        forged = copy.deepcopy(curve)
        forged['rows'][0]['cumulative']['completed']['numerator'] = 1
        with self.assertRaises(ValueError):
            verify_curve(forged)

    def test_naive_and_aligned_estimates_expose_the_bias(self):
        observations = [episode('a', 1, 'completed'), episode('b', 2, 'death'),
                        episode('c', 3, 'transfer'), episode('d', 4, None)]
        curve = cumulative_incidence(observations)
        aligned = aligned_naive(observations, curve)
        self.assertEqual(len(aligned['rows']), len(curve['rows']))
        self.assertFalse(aligned['for_conclusion'])
        self.assertEqual(aligned['rows'][0]['aalen_johansen_completed'], curve['rows'][0]['cumulative']['completed'])
        last = aligned['rows'][-1]
        self.assertGreaterEqual(as_fraction(last['naive_completed']),
                                as_fraction(last['aalen_johansen_completed']))
        self.assertEqual(as_fraction(last['naive_completed']), Fraction(1, 4))
        self.assertEqual(as_fraction(last['aalen_johansen_completed']), Fraction(1, 4))
        naive = naive_estimates(observations, 2)
        self.assertEqual(naive['targets']['completed_as_event_other_discharges_censored']['cumulative'],
                         {'numerator': 1, 'denominator': 4})
        self.assertEqual(naive['targets']['any_discharge_as_event']['cumulative'],
                         {'numerator': 1, 'denominator': 2})
        self.assertFalse(naive['for_conclusion'])
        self.assertEqual(naive['formal_estimand'], 'aalen_johansen_cumulative_incidence')

    def test_target_beyond_followup_is_not_extrapolated(self):
        curve = cumulative_incidence([episode('a', 2, 'completed'), episode('b', 5, None)])
        self.assertFalse(target_at(curve, 6)['estimable'])
        self.assertIsNone(target_at(curve, 6)['survival'])
        landed = target_at(cumulative_incidence([episode('a', 2, 'completed')]), 40)
        self.assertTrue(landed['estimable'])
        self.assertEqual(as_fraction(landed['survival']), 0)
        empty = cumulative_incidence([])
        self.assertFalse(target_at(empty, 1)['estimable'])
        self.assertEqual(target_at(empty, 1)['reason'], 'no_observations')



    def test_aligned_naive_diverges_when_competing_events_are_censored_away(self):
        observations = [episode('a', 1, 'death'), episode('b', 2, 'death'),
                        episode('c', 3, 'completed'), episode('d', 4, 'transfer')]
        design = micro([
            visit('b1', 's1', '2026-11-01', '2026-11-02', 'death'),
            visit('b2', 's2', '2026-11-01', '2026-11-03', 'death'),
            visit('b3', 's3', '2026-11-01', '2026-11-04', 'completed'),
            visit('b4', 's4', '2026-11-01', '2026-11-05', 'transfer')],
            people=('s1', 's2', 's3', 's4'))
        dataset = Dataset(design)
        report = analyze(dataset, '2026-11-12', '2026-11-12', study_stage='acute',
                         study_target_day=4)
        study = report['competing_risks']
        self.assertEqual(as_fraction(study['target']['cumulative']['completed']), Fraction(1, 4))
        aligned = study['aligned_naive']['rows']
        self.assertEqual(len(aligned), len(study['curve']['rows']))
        self.assertEqual(aligned[-1]['time'], 4)
        last = aligned[-1]
        self.assertEqual(as_fraction(last['aalen_johansen_completed']), Fraction(1, 4))
        self.assertGreater(as_fraction(last['naive_completed']), Fraction(1, 4))
        self.assertEqual(as_fraction(last['naive_completed']), Fraction(1, 2))
        naive = study['naive']['targets']['any_discharge_as_event']
        self.assertTrue(naive['estimable'])
        self.assertEqual(as_fraction(naive['cumulative']), 1)
        self.assertFalse(study['naive']['for_conclusion'])

class CompetingRiskSelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.data = micro([
            visit('pre1', 'q1', '2026-10-30', '2026-11-01', 'transfer'),
            visit('pre2', 'q1', '2026-11-03', '2026-11-05', 'completed'),
            visit('open1', 'q2', '2026-10-31', None, None),
            visit('new1', 'q3', '2026-11-03', '2026-11-06', 'death'),
            visit('new2', 'q3', '2026-11-07', '2026-11-07', 'completed'),
            visit('open2', 'q4', '2026-11-04', None, None),
            visit('late1', 'q5', '2026-11-04', '2026-11-05', 'completed', recorded='2026-11-13'),
            visit('pre3', 'q6', '2026-11-02', '2026-11-02', 'left')],
            people=('q1', 'q2', 'q3', 'q4', 'q5', 'q6'))
        self.dataset = Dataset(self.data)
        self.report = analyze(self.dataset, '2026-11-12', '2026-11-12', horizon=7,
                              study_stage='acute', study_window_start='2026-11-03',
                              study_target_day=7)
        self.study = self.report['competing_risks']

    def tearDown(self):
        self.temp.cleanup()

    def test_first_visible_episode_per_person_and_exclusion_reasons(self):
        selected = self.study['selected']
        self.assertEqual([row['visit'] for row in self.study['episodes']], ['new1', 'open2'])
        self.assertEqual(selected['window']['start'], '2026-11-03')
        reasons = {row['visit']: row['reason'] for row in selected['excluded_records']}
        self.assertEqual(reasons, {'pre1': 'admitted_before_window_start',
                                   'open1': 'admitted_before_window_start',
                                   'pre2': 'after_prior_stage_admission_before_window_start',
                                   'new2': 'repeat_admission_after_first_episode',
                                   'late1': 'recorded_after_knowledge_cutoff',
                                   'pre3': 'admitted_before_window_start'})
        self.assertEqual(selected['exclusion_counts']['records'],
                         {'admitted_before_window_start': 3,
                          'after_prior_stage_admission_before_window_start': 1,
                          'recorded_after_knowledge_cutoff': 1,
                          'repeat_admission_after_first_episode': 1})
        self.assertEqual(selected['exclusion_counts']['people'],
                         {'already_inpatient_at_window_start': 1,
                          'no_eligible_visible_admission_in_window': 1,
                          'prior_stage_admission_before_window_start': 2})
        self.assertEqual([row['visit'] for row in selected['invisible_records']], ['late1'])
        self.assertTrue(verify_study(self.dataset, self.dataset.slice('2026-11-12', '2026-11-12'),
                                     'acute', self.study)['valid'])

    def test_censoring_and_curve_arithmetic(self):
        curve = self.study['curve']
        self.assertEqual(curve['population'], 2)
        self.assertEqual(curve['censored'], 1)
        self.assertEqual(curve['events'], {'completed': 0, 'transfer': 0, 'death': 1, 'left': 0})
        self.assertEqual(as_fraction(curve['rows'][0]['cumulative']['death']), Fraction(1, 2))
        self.assertEqual(curve['rows'][-1]['time'], 8)
        target = self.study['target']
        self.assertTrue(target['estimable'])
        self.assertEqual(as_fraction(target['cumulative']['death']), Fraction(1, 2))
        self.assertEqual(as_fraction(target['cumulative']['completed']), 0)
        self.assertEqual(as_fraction(target['survival']), Fraction(1, 2))
        self.assertTrue(verify(self.dataset, self.report, True)['replayed'])

    def test_recording_later_creates_another_knowledge_version(self):
        late = micro([
            visit('k1', 'r1', '2026-11-04', '2026-11-05', 'completed'),
            visit('k2', 'r2', '2026-11-04', '2026-11-06', 'death'),
            visit('k3', 'r3', '2026-11-04', '2026-11-05', 'completed', recorded='2026-11-13')],
            people=('r1', 'r2', 'r3'))
        late_dataset = Dataset(late)
        early_study = analyze(late_dataset, '2026-11-12', '2026-11-12', study_stage='acute',
                              study_target_day=7)['competing_risks']
        later_study = analyze(late_dataset, '2026-11-12', '2026-11-13', study_stage='acute',
                              study_target_day=7)['competing_risks']
        self.assertEqual(early_study['curve']['population'], 2)
        self.assertEqual(later_study['curve']['population'], 3)
        self.assertIn('k3', [row['visit'] for row in early_study['selected']['invisible_records']])
        self.assertNotIn('k3', [row['visit'] for row in later_study['selected']['invisible_records']])
        self.assertEqual(as_fraction(early_study['target']['cumulative']['completed']), Fraction(1, 2))
        self.assertEqual(as_fraction(later_study['target']['cumulative']['completed']), Fraction(2, 3))
        self.assertTrue(verify_study(late_dataset, late_dataset.slice('2026-11-12', '2026-11-13'),
                                     'acute', later_study)['valid'])
        with self.assertRaises(ValueError):
            verify_study(late_dataset, late_dataset.slice('2026-11-12', '2026-11-12'),
                         'acute', later_study)

    def test_future_business_discharge_is_censored_not_counted(self):
        future = analyze(self.dataset, '2026-11-05', '2026-11-12', study_stage='acute',
                         study_window_start='2026-11-03', study_target_day=2)['competing_risks']
        episode = next(row for row in future['episodes'] if row['visit'] == 'new1')
        self.assertIsNone(episode['outcome'])
        self.assertEqual(episode['duration'], 2)
        self.assertEqual(episode['censoring_reason'], 'discharge_after_business_cutoff')
        self.assertEqual(as_fraction(future['target']['cumulative']['death']), 0)
        self.assertEqual(future['curve']['events']['death'], 0)
        self.assertTrue(all(row['event'] is False for row in future['episodes']))
        self.assertEqual(future['curve']['population'], 2)


class CompetingRiskComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.data = micro([
            visit('a1', 'e1', '2026-11-01', '2026-11-02', 'completed'),
            visit('a2', 'e2', '2026-11-01', '2026-11-03', 'death'),
            visit('a3', 'e3', '2026-11-01', None, None),
            visit('b1', 'w1', '2026-11-01', '2026-11-02', 'transfer'),
            visit('b2', 'w2', '2026-11-01', '2026-11-04', 'completed'),
            visit('b3', 'w3', '2026-11-01', '2026-11-04', 'completed')],
            people=('e1', 'e2', 'e3', 'w1', 'w2', 'w3'))
        self.data['people'][0].update(county='east')
        for index, person in enumerate(('e1', 'e2', 'e3')):
            self.data['people'][index]['county'] = 'east'
        for index, person in enumerate(('w1', 'w2', 'w3'), start=3):
            self.data['people'][index]['county'] = 'west'
        self.dataset = Dataset(self.data)
        self.spec = {'id': 'county-completed', 'dimension': 'county', 'left': ['east'],
                     'right': ['west'], 'outcome': 'completed', 'replicates': 60, 'seed': 11}

    def tearDown(self):
        self.temp.cleanup()

    def test_risk_difference_and_reproducible_person_bootstrap(self):
        report = analyze(self.dataset, '2026-11-12', '2026-11-12', study_stage='acute',
                         study_target_day=4, study_comparisons=[self.spec])
        row = report['competing_risks']['comparisons'][0]
        self.assertEqual(as_fraction(row['left']['cumulative']['completed']), Fraction(1, 3))
        self.assertEqual(as_fraction(row['right']['cumulative']['completed']), Fraction(2, 3))
        self.assertEqual(as_fraction(row['risk_difference']), Fraction(-1, 3))
        bootstrap = row['bootstrap']
        self.assertEqual(bootstrap['replicates'], 60)
        self.assertEqual(bootstrap['supported_replicates'] + bootstrap['unsupported_replicates'], 60)
        self.assertTrue(bootstrap['unsupported_replicates_retained'])
        self.assertLessEqual(bootstrap['interval'][0], bootstrap['interval'][1])
        replay = analyze(self.dataset, '2026-11-12', '2026-11-12', study_stage='acute',
                         study_target_day=4, study_comparisons=[self.spec])
        self.assertEqual(replay, report)

    def test_unsupported_target_days_remain_counted_in_bootstrap(self):
        spec = dict(self.spec, outcome='death')
        report = analyze(self.dataset, '2026-11-12', '2026-11-12', study_stage='acute',
                         study_target_day=12, study_comparisons=[spec])
        row = report['competing_risks']['comparisons'][0]
        self.assertFalse(row['left']['estimable'])
        self.assertIsNone(row['risk_difference'])
        bootstrap = row['bootstrap']
        self.assertEqual(bootstrap['supported_replicates'] + bootstrap['unsupported_replicates'], 60)
        self.assertGreater(bootstrap['unsupported_replicates'], 0)
        self.assertTrue(bootstrap['unsupported_replicates_retained'])
        self.assertLessEqual(bootstrap['supported_replicates'], 60)

    def test_fixed_cohorts_and_overlap_are_rejected(self):
        with self.assertRaises(ValueError):
            analyze(self.dataset, '2026-11-12', '2026-11-12', study_stage='acute',
                    study_comparisons=[dict(self.spec, left=['east'], right=['east'])])
        with self.assertRaises(ValueError):
            analyze(self.dataset, '2026-11-12', '2026-11-12', study_stage='acute',
                    study_comparisons=[dict(self.spec, id='a'), dict(self.spec, id='a')])
        with self.assertRaises(ValueError):
            analyze(self.dataset, '2026-11-12', '2026-11-12', study_stage='acute',
                    study_comparisons=[dict(self.spec, outcome='recovered')])
        with self.assertRaises(ValueError):
            analyze(self.dataset, '2026-11-12', '2026-11-12', study_stage='ward')
        narrow = analyze(self.dataset, '2026-11-12', '2026-11-12', study_stage='acute',
                         study_comparisons=[self.spec])
        self.assertEqual(narrow['competing_risks']['definition']['comparisons'], [self.spec])

    def test_report_checksum_rejects_tampered_competing_curve(self):
        report = analyze(self.dataset, '2026-11-12', '2026-11-12', study_stage='acute',
                         study_target_day=4, study_comparisons=[self.spec])
        forged = copy.deepcopy(report)
        forged['competing_risks']['curve']['rows'][0]['cumulative']['death']['numerator'] += 1
        with self.assertRaises(ValueError):
            verify(self.dataset, forged)
        self.assertTrue(verify(self.dataset, report, True)['valid'])


class CompetingRiskIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.data = json.loads(Path('fixtures/care.json').read_text())
        self.dataset = Dataset(self.data)
        self.options = {'as_of': '2026-11-12', 'known_at': '2026-11-12', 'horizon': 7,
                        'study_stage': 'acute', 'study_target_day': 6}
        self.report = analyze(self.dataset, **self.options)

    def tearDown(self):
        self.temp.cleanup()

    def test_fixture_study_and_plan_dispatch(self):
        study = self.report['competing_risks']
        self.assertEqual([row['visit'] for row in study['episodes']], ['v1', 'v5', 'v4', 'v7'])
        self.assertEqual(as_fraction(study['target']['cumulative']['death']), Fraction(1, 4))
        self.assertEqual(as_fraction(study['target']['cumulative']['transfer']), Fraction(1, 2))
        self.assertEqual(as_fraction(study['target']['cumulative']['completed']), Fraction(1, 4))
        definition = {'version': 1, 'operations': [
            {'id': 'curve', 'kind': 'competing-risks', 'options': {}},
            {'id': 'day', 'kind': 'competing-target', 'options': {'day': 6}},
            {'id': 'visit', 'kind': 'competing-lineage', 'options': {'visit': 'v1'}}]}
        execution = run_plan(self.dataset, self.report, definition)
        self.assertTrue(verify_plan(self.dataset, self.report, execution)['valid'])
        results = {row['id']: row['result'] for row in execution['results']}
        self.assertEqual(results['day'], study['target'])
        self.assertTrue(results['visit']['selected'])
        self.assertEqual(results['visit']['episode']['outcome'], 'transfer')
        self.assertTrue(results['visit']['term']['cumulative_incidence_terms_for_episode'])
        with self.assertRaises(ValueError):
            run_plan(self.dataset, self.report, {'version': 1, 'operations': [
                {'id': 'bad', 'kind': 'competing-compare', 'options': {'comparison': 'absent'}}]})

    def test_query_families_and_lineage_for_excluded_record(self):
        self.assertEqual(dispatch(self.dataset, self.report, 'competing-risks'),
                         self.report['competing_risks'])
        target = dispatch(self.dataset, self.report, 'competing', {'target_day': 1})
        self.assertTrue(target['estimable'])
        with self.assertRaises(ValueError):
            dispatch(self.dataset, self.report, 'competing', {'target_day': 1, 'comparison': 'none'})
        lineage = competing_lineage(self.dataset, self.report, 'v6')
        self.assertFalse(lineage['selected'])
        self.assertTrue(lineage['outside_selected_stage'])
        self.assertIsNone(lineage['excluded'])
        self.assertFalse(lineage['counts_to_curve'])
        self.assertIn('records', lineage)
        with self.assertRaises(KeyError):
            competing_lineage(self.dataset, self.report, 'absent')

    def test_journal_archive_and_restatement_keep_knowledge_versions(self):
        with Journal(self.root / 'care.db') as journal:
            journal.replace('first', 0, self.data)
            record = journal.analyze('competing', 1, self.options)
            self.assertEqual(record['result']['competing_risks'], self.report['competing_risks'])
            self.assertTrue(journal.audit()['valid'])
            options = dict(self.options, study_stage='primary', study_target_day=2)
            first = journal.analyze('primary-open', 1, options)
            self.assertEqual(first['result']['competing_risks']['curve']['censored'], 1)
            changed = copy.deepcopy(self.data)
            next(row for row in changed['visits'] if row['id'] == 'v8')['discharged'] = '2026-11-12'
            next(row for row in changed['visits'] if row['id'] == 'v8')['outcome'] = 'completed'
            next(row for row in changed['visits'] if row['id'] == 'v8')['recorded'] = '2026-11-12'
            journal.replace('closed', 1, changed)
            later = journal.analyze('primary-closed', 2, options)
            self.assertNotEqual(later['result']['competing_risks']['curve'],
                                first['result']['competing_risks']['curve'])
            self.assertEqual(journal.analysis('primary-open'), first)

    def test_cli_options_are_part_of_analysis_identity(self):
        with Journal(self.root / 'care.db') as journal:
            journal.replace('first', 0, self.data)
            first = journal.analyze('one', 1, self.options)
            from encounters.values import Conflict
            with self.assertRaises(Conflict):
                journal.analyze('one', 1, dict(self.options, study_target_day=9))
            dispatch_result = dispatch(self.dataset, first['result'], 'competing', {'target_day': 2})
            self.assertTrue(dispatch_result['estimable'])


if __name__ == '__main__':
    unittest.main()

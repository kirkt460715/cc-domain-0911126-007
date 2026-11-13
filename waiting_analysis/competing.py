"""Aalen-Johansen competing-risk study of inpatient discharge destinations.

One visible inpatient episode per person and stage is selected.  Time zero is
admission; an episode still open at the business cutoff is right-censored.
The four discharge destinations are mutually exclusive.  Estimates are
descriptive and do not support clinical decisions.
"""
from collections import defaultdict
from fractions import Fraction
import random
from encounters.values import day, days, digest, fields, identifier, integer, ratio, Invalid
from encounters.schema import STAGES
from waiting_analysis.survival import kaplan_meier, survival_at


OUTCOMES = ('completed', 'transfer', 'death', 'left')
DIMENSIONS = ('county', 'age_band')


def fraction(value):
    return {'numerator': value.numerator, 'denominator': value.denominator}


def as_fraction(record):
    if not isinstance(record, dict) or set(record) != {'numerator', 'denominator'}:
        raise Invalid('fraction record required')
    integer(record['numerator'], 'numerator', -1000000000, 1000000000)
    integer(record['denominator'], 'denominator', 1, 1000000000)
    return Fraction(record['numerator'], record['denominator'])


def validate_episodes(episodes):
    if not isinstance(episodes, list) or len(episodes) > 100000:
        raise Invalid('bounded episodes required')
    keys = set()
    for row in episodes:
        fields(row, ('id', 'duration', 'outcome'))
        if not isinstance(row['id'], str) or not row['id'] or row['id'] in keys:
            raise Invalid('duplicate or invalid episode identity')
        keys.add(row['id'])
        integer(row['duration'], 'duration', 0, 36500)
        if row['outcome'] is not None and row['outcome'] not in OUTCOMES:
            raise Invalid('unsupported discharge outcome')
    return episodes


def cumulative_incidence(episodes):
    """Aalen-Johansen cumulative incidence; same-day events precede censoring."""
    validate_episodes(episodes)
    grouped = {}
    totals = {key: 0 for key in OUTCOMES}
    censored_total = 0
    for row in episodes:
        counts = grouped.setdefault(row['duration'], {'outcomes': {key: 0 for key in OUTCOMES}, 'censored': 0})
        if row['outcome'] is None:
            counts['censored'] += 1
            censored_total += 1
        else:
            counts['outcomes'][row['outcome']] += 1
            totals[row['outcome']] += 1
    at_risk = len(episodes)
    survival = Fraction(1)
    cumulative = {key: Fraction(0) for key in OUTCOMES}
    rows = []
    for time, counts in sorted(grouped.items()):
        events = sum(counts['outcomes'].values())
        censored = counts['censored']
        if events > at_risk or events + censored > at_risk:
            raise Invalid('event and censoring counts exceed risk set')
        before = survival
        if events:
            survival *= Fraction(at_risk - events, at_risk)
            for key in OUTCOMES:
                if counts['outcomes'][key]:
                    cumulative[key] += before * Fraction(counts['outcomes'][key], at_risk)
        total_cumulative = sum(cumulative.values(), Fraction(0))
        if survival + total_cumulative != 1:
            raise Invalid('competing-risk probability conservation failed')
        rows.append({'time': time, 'at_risk': at_risk, 'event_time': bool(events),
                     'outcome_events': dict(counts['outcomes']), 'events_total': events,
                     'censored': censored, 'survival_before': fraction(before),
                     'survival': fraction(survival),
                     'cumulative': {key: fraction(cumulative[key]) for key in OUTCOMES},
                     'cumulative_total': fraction(total_cumulative),
                     'conservation': fraction(survival + total_cumulative)})
        at_risk -= events + censored
    if at_risk != 0:
        raise Invalid('risk sets do not conserve episodes')
    return {'population': len(episodes), 'outcomes': list(OUTCOMES),
            'events': dict(totals), 'events_total': sum(totals.values()),
            'censored': censored_total, 'rows': rows,
            'last_time': rows[-1]['time'] if rows else None,
            'last_survival': fraction(survival),
            'ties': 'events_before_censoring',
            'definition': 'aalen_johansen_cumulative_incidence',
            'rows_include_censoring_only_times': True}


def target_at(curve, time):
    """Estimate at a target day; refuse to project past supported follow-up.

    Support ends at the last observed time unless every episode had already
    reached an outcome by then (survival zero), so a later day is estimable.
    """
    integer(time, 'target_day', 0, 36500)
    if curve['population'] == 0 or not curve['rows']:
        return {'day': time, 'estimable': False, 'reason': 'no_observations',
                'survival': None, 'cumulative': None}
    last = curve['rows'][-1]
    supported = time <= last['time'] or as_fraction(last['survival']) == 0
    if not supported:
        return {'day': time, 'estimable': False, 'reason': 'followup_ends_before_target',
                'survival': None, 'cumulative': None,
                'last_supported_day': last['time']}
    survival = Fraction(1)
    cumulative = {key: Fraction(0) for key in OUTCOMES}
    for row in curve['rows']:
        if row['time'] > time:
            break
        survival = as_fraction(row['survival'])
        cumulative = {key: as_fraction(row['cumulative'][key]) for key in OUTCOMES}
    return {'day': time, 'estimable': True, 'reason': None,
            'survival': fraction(survival),
            'cumulative': {key: fraction(cumulative[key]) for key in OUTCOMES},
            'cumulative_total': fraction(survival + sum(cumulative.values(), Fraction(0))),
            'last_supported_day': last['time']}


def verify(curve):
    if not isinstance(curve, dict) or curve.get('outcomes') != list(OUTCOMES):
        raise Invalid('unsupported competing-risk curve')
    at_risk = curve['population']
    previous_time = -1
    previous_survival = Fraction(1)
    cumulative = {key: Fraction(0) for key in OUTCOMES}
    events = {key: 0 for key in OUTCOMES}
    censored = 0
    for row in curve['rows']:
        if row['time'] <= previous_time or row['at_risk'] != at_risk:
            raise Invalid('competing-risk curve ordering or risk set invalid')
        counts = row['outcome_events']
        if set(counts) != set(OUTCOMES) or any(type(counts[key]) is not int or counts[key] < 0 for key in OUTCOMES):
            raise Invalid('competing-risk outcome counts invalid')
        total_events = sum(counts.values())
        if total_events != row['events_total'] or type(row['censored']) is not int or row['censored'] < 0:
            raise Invalid('competing-risk event total invalid')
        if as_fraction(row['survival_before']) != previous_survival:
            raise Invalid('competing-risk survival before event mismatch')
        expected = previous_survival * Fraction(at_risk - total_events, at_risk)
        if as_fraction(row['survival']) != expected:
            raise Invalid('competing-risk survival arithmetic mismatch')
        for key in OUTCOMES:
            cumulative[key] += previous_survival * Fraction(counts[key], at_risk)
            if as_fraction(row['cumulative'][key]) != cumulative[key]:
                raise Invalid('competing-risk cumulative incidence mismatch')
        total = sum(cumulative.values(), Fraction(0))
        if as_fraction(row['cumulative_total']) != total or as_fraction(row['conservation']) != expected + total:
            raise Invalid('competing-risk conservation record mismatch')
        if expected + total != 1:
            raise Invalid('competing-risk probabilities do not sum to one')
        at_risk -= total_events + row['censored']
        for key in OUTCOMES:
            events[key] += counts[key]
        censored += row['censored']
        previous_time = row['time']
        previous_survival = expected
    if at_risk != 0:
        raise Invalid('competing-risk risk sets do not conserve episodes')
    if events != curve['events'] or censored != curve['censored']:
        raise Invalid('competing-risk curve totals mismatch')
    if curve['rows'] and previous_survival != as_fraction(curve['last_survival']):
        raise Invalid('competing-risk final survival mismatch')
    return {'valid': True, 'population': curve['population'], 'rows': len(curve['rows'])}


def naive_estimates(episodes, target_day):
    """Kaplan-Meier estimates that deliberately ignore the competing structure."""
    completion = [{'id': row['id'], 'duration': row['duration'],
                   'event': row['outcome'] == 'completed'} for row in episodes]
    any_discharge = [{'id': row['id'], 'duration': row['duration'],
                      'event': row['outcome'] is not None} for row in episodes]
    curves = {'completed_as_event_other_discharges_censored': kaplan_meier(completion),
              'any_discharge_as_event': kaplan_meier(any_discharge)}
    targets = {}
    for name, curve in curves.items():
        if curve['population'] == 0 or not curve['rows']:
            targets[name] = {'day': target_day, 'estimable': False, 'reason': 'no_observations',
                             'cumulative': None}
            continue
        last = curve['rows'][-1]
        if target_day > last['time'] and as_fraction(last['survival']) != 0:
            targets[name] = {'day': target_day, 'estimable': False,
                             'reason': 'followup_ends_before_target',
                             'last_supported_day': last['time'], 'cumulative': None}
        else:
            targets[name] = {'day': target_day, 'estimable': True, 'reason': None,
                             'cumulative': fraction(1 - survival_at(curve, target_day))}
    return {'curves': curves, 'targets': targets,
            'formal_estimand': 'aalen_johansen_cumulative_incidence',
            'for_conclusion': False,
            'warning': 'naive_kaplan_meier_censors_competing_discharges_and_then_uses_one_minus_survival'}


def aligned_naive(episodes, curve):
    """Place the naive one-minus-KM completion estimate beside every event time."""
    observations = [{'id': row['id'], 'duration': row['duration'],
                     'event': row['outcome'] == 'completed'} for row in episodes]
    naive = kaplan_meier(observations)
    rows = []
    for row in curve['rows']:
        time = row['time']
        if naive['population'] == 0 or not naive['rows']:
            value = None
        else:
            last = naive['rows'][-1]
            supported = time <= last['time'] or as_fraction(last['survival']) == 0
            value = fraction(1 - survival_at(naive, time)) if supported else None
        rows.append({'time': time,
                     'aalen_johansen_completed': row['cumulative']['completed'],
                     'naive_completed': value})
    return {'rows': rows, 'definition': 'naive_columns_ignore_competing_discharges',
            'for_conclusion': False}


def _visible_and_invisible(dataset, view, stage):
    business, knowledge = day(view['as_of']), day(view['known_at'])
    visible, invisible = [], []
    for key, stored in sorted(dataset.visits.items()):
        if dataset.stage(key) != stage:
            continue
        if day(stored['recorded']) > knowledge:
            invisible.append({'visit': key, 'person': stored['person'],
                              'reason': 'recorded_after_knowledge_cutoff'})
        elif day(stored['admitted']) > business:
            invisible.append({'visit': key, 'person': stored['person'],
                              'reason': 'admitted_after_business_cutoff'})
        else:
            visible.append((key, view['visits'][key]))
    return visible, invisible


def _implied_start(visible, business):
    """Start of the observation window covering every visible stage admission."""
    admissions = [day(row['admitted']) for _, row in visible]
    return min(admissions, default=business).isoformat()


def selection(dataset, view, stage, window_start=None):
    """Select the first visible episode per person in one stage's window."""
    if stage not in STAGES:
        raise Invalid('unsupported stage filter')
    business = day(view['as_of'])
    visible, invisible = _visible_and_invisible(dataset, view, stage)
    start = window_start
    if start is None:
        start = _implied_start(visible, business)
    first = day(start, 'window_start')
    if first > business:
        raise Invalid('observation window starts after business cutoff')
    visible_by_person = defaultdict(list)
    for key, row in visible:
        visible_by_person[row['person']].append((key, row))
    invisible_by_person = defaultdict(list)
    for row in invisible:
        invisible_by_person[row['person']].append(row)
    people = sorted(set(visible_by_person) | set(invisible_by_person))
    episodes, excluded_records, excluded_people = [], [], []
    for person in people:
        person_visible = sorted(visible_by_person[person], key=lambda item: (item[1]['admitted'], item[0]))
        before = [(key, row) for key, row in person_visible if day(row['admitted']) < first]
        candidates = [(key, row) for key, row in person_visible if day(row['admitted']) >= first]
        for key, row in before:
            excluded_records.append({'visit': key, 'person': person,
                                     'reason': 'admitted_before_window_start',
                                     'admitted': row['admitted'], 'discharged': row['discharged'],
                                     'recorded': row['recorded']})
        for row in invisible_by_person[person]:
            excluded_records.append({'visit': row['visit'], 'person': person,
                                     'reason': row['reason']})
        if before:
            ongoing = [key for key, row in before
                       if row['discharged'] is None or day(row['discharged']) > first]
            reason = 'already_inpatient_at_window_start' if ongoing else 'prior_stage_admission_before_window_start'
            record_reason = ('person_already_inpatient_at_window_start' if ongoing
                             else 'after_prior_stage_admission_before_window_start')
            excluded_people.append({'person': person, 'reason': reason,
                                    'prewindow_visits': [key for key, _ in before],
                                    'ongoing_visits': ongoing,
                                    'records': [key for key, _ in candidates]})
            for key, row in candidates:
                excluded_records.append({'visit': key, 'person': person,
                                         'reason': record_reason,
                                         'admitted': row['admitted'], 'discharged': row['discharged'],
                                         'recorded': row['recorded']})
            continue
        if not candidates:
            excluded_people.append({'person': person, 'reason': 'no_eligible_visible_admission_in_window',
                                    'records': [key for key, item in person_visible
                                                if day(item['admitted']) < first]
                                    + [row['visit'] for row in invisible_by_person[person]]})
            continue
        chosen, row = candidates[0]
        source = view['visits'][chosen]
        discharged = source['discharged']
        observed = source['outcome'] if discharged is not None else None
        censoring_reason = None
        if discharged is None:
            stored = dataset.visits[chosen]
            censoring_reason = ('discharge_after_business_cutoff'
                                if stored['discharged'] and day(stored['discharged']) > business
                                else 'open_at_business_cutoff')
        episodes.append({'visit': chosen, 'person': person, 'stage': stage,
                         'facility': source['facility'],
                         'county': dataset.people[person]['county'],
                         'age_band': dataset.people[person]['age_band'],
                         'rural': dataset.people[person]['rural'],
                         'admitted': source['admitted'], 'discharged': discharged,
                         'outcome': observed,
                         'duration': days(source['admitted'], discharged or view['as_of']),
                         'event': observed is not None,
                         'recorded': source['recorded'],
                         'source_digest': digest(source),
                         'censoring_reason': censoring_reason,
                         'first_visible_episode': True})
        for key, row in candidates[1:]:
            excluded_records.append({'visit': key, 'person': person,
                                     'reason': 'repeat_admission_after_first_episode',
                                     'first_visit': chosen, 'admitted': row['admitted'],
                                     'discharged': row['discharged'], 'recorded': row['recorded']})
    episodes.sort(key=lambda row: (row['admitted'], row['visit']))
    record_counts = defaultdict(int)
    for row in excluded_records:
        record_counts[row['reason']] += 1
    person_counts = defaultdict(int)
    for row in excluded_people:
        person_counts[row['reason']] += 1
    return {'stage': stage,
            'window': {'start': start, 'end': view['as_of'],
                       'rule': 'admission_on_or_after_start_and_on_or_before_business_cutoff'},
            'time_origin': 'admission',
            'episodes': episodes,
            'invisible_records': invisible,
            'excluded_records': excluded_records,
            'excluded_people': excluded_people,
            'exclusion_counts': {'records': dict(sorted(record_counts.items())),
                                 'people': dict(sorted(person_counts.items()))},
            'stage_records': len(visible) + len(invisible),
            'visible_stage_records': len(visible),
            'one_episode_per_person': True,
            'first_episode_rule': 'earliest_visible_admission_per_person_within_selected_stage',
            'repeat_admission_rule': 'later_visible_admissions_in_same_stage_excluded',
            'prewindow_rule': 'people_with_a_stage_admission_before_window_start_excluded_as_not_new_admissions',
            'ongoing_rule': 'people_inpatient_at_window_start_reported_separately',
            'non_visible_rule': 'records_registered_after_known_at_never_enter_the_cohort_and_do_not_shape_the_window'}


def _subject(episodes):
    return [{'id': row['visit'], 'duration': row['duration'], 'outcome': row['outcome']}
            for row in episodes]


def _bootstrap_comparison(left, right, outcome, target_day, replicates, seed):
    if not left or not right:
        return {'replicates': replicates, 'seed': seed, 'supported_replicates': 0,
                'unsupported_replicates': replicates, 'interval': None,
                'reason': 'empty_cohort', 'unsupported_replicates_retained': True}
    generator = random.Random(seed)
    values = []
    unsupported = 0
    for index in range(replicates):
        sampled = []
        for side, rows in (('l', left), ('r', right)):
            sampled.append([{'id': str(index) + ':' + side + ':' + str(draw),
                             'duration': original['duration'], 'outcome': original['outcome']}
                            for draw, original in enumerate(generator.choices(rows, k=len(rows)))])
        left_target = target_at(cumulative_incidence(sampled[0]), target_day)
        right_target = target_at(cumulative_incidence(sampled[1]), target_day)
        if not left_target['estimable'] or not right_target['estimable']:
            unsupported += 1
            continue
        values.append(float(as_fraction(left_target['cumulative'][outcome])
                            - as_fraction(right_target['cumulative'][outcome])))
    values.sort()
    interval = None
    if values:
        interval = [values[int((len(values) - 1) * .025)], values[int((len(values) - 1) * .975)]]
    return {'replicates': replicates, 'seed': seed,
            'supported_replicates': len(values), 'unsupported_replicates': unsupported,
            'unsupported_replicates_retained': True,
            'supported_fraction': ratio(len(values), replicates),
            'interval': interval,
            'interval_method': 'percentile_over_supported_replicates_only',
            'warning': 'unsupported_replicates_are_counted_not_dropped'}


def comparison(episodes, spec, target_day):
    fields(spec, ('id', 'dimension', 'left', 'right', 'outcome', 'replicates', 'seed'))
    identifier(spec['id'], 'comparison.id')
    if spec['dimension'] not in DIMENSIONS:
        raise Invalid('unsupported comparison dimension')
    if spec['outcome'] not in OUTCOMES:
        raise Invalid('unsupported comparison outcome')
    for key in ('left', 'right'):
        values = spec[key]
        if not isinstance(values, list) or not values or len(values) > 100:
            raise Invalid('comparison cohort must be a bounded nonempty array')
        if any(not isinstance(value, str) for value in values):
            raise Invalid('comparison cohort value must be text')
    if set(spec['left']) & set(spec['right']):
        raise Invalid('comparison cohorts overlap by definition')
    integer(spec['replicates'], 'replicates', 20, 2000)
    integer(spec['seed'], 'seed', 0, 2 ** 32 - 1)
    left = [row for row in episodes if row[spec['dimension']] in spec['left']]
    right = [row for row in episodes if row[spec['dimension']] in spec['right']]
    if {row['person'] for row in left} & {row['person'] for row in right}:
        raise Invalid('comparison cohorts overlap by person')
    left_curve = cumulative_incidence(_subject(left))
    right_curve = cumulative_incidence(_subject(right))
    left_target = target_at(left_curve, target_day)
    right_target = target_at(right_curve, target_day)
    left_aligned = aligned_naive(_subject(left), left_curve)
    right_aligned = aligned_naive(_subject(right), right_curve)
    difference = None
    if left_target['estimable'] and right_target['estimable']:
        difference = fraction(as_fraction(left_target['cumulative'][spec['outcome']])
                              - as_fraction(right_target['cumulative'][spec['outcome']]))
    def naive(rows):
        observations = [{'id': row['visit'], 'duration': row['duration'],
                         'event': row['outcome'] == spec['outcome']} for row in rows]
        curve = kaplan_meier(observations)
        if curve['population'] == 0 or not curve['rows']:
            return None
        last = curve['rows'][-1]
        if target_day > last['time'] and as_fraction(last['survival']) != 0:
            return None
        return fraction(1 - survival_at(curve, target_day))
    left_naive, right_naive = naive(left), naive(right)
    naive_difference = (fraction(as_fraction(left_naive) - as_fraction(right_naive))
                        if left_naive is not None and right_naive is not None else None)
    return {'id': spec['id'], 'dimension': spec['dimension'],
            'left_cohort': spec['left'], 'right_cohort': spec['right'],
            'outcome': spec['outcome'], 'target_day': target_day,
            'left': {'people': len({row['person'] for row in left}), 'episodes': len(left),
                     'estimable': left_target['estimable'], 'reason': left_target['reason'],
                     'cumulative': left_target['cumulative'],
                     'naive_cumulative': left_naive,
                     'aligned_naive': left_aligned},
            'right': {'people': len({row['person'] for row in right}), 'episodes': len(right),
                      'estimable': right_target['estimable'], 'reason': right_target['reason'],
                      'cumulative': right_target['cumulative'],
                      'naive_cumulative': right_naive,
                      'aligned_naive': right_aligned},
            'risk_difference': difference,
            'naive_risk_difference': naive_difference,
            'bootstrap': _bootstrap_comparison(left, right, spec['outcome'],
                                               target_day, spec['replicates'], spec['seed']),
            'formal_estimand': 'aalen_johansen_cumulative_incidence',
            'causal_interpretation': False,
            'cohorts_fixed_by_study': True}


def study(dataset, view, stage, window_start=None, target_day=None, comparisons=None):
    if stage not in STAGES:
        raise Invalid('unsupported stage filter')
    target = view['as_of'] if target_day is None else target_day
    integer(target, 'target_day', 0, 36500)
    comparisons = [] if comparisons is None else comparisons
    if not isinstance(comparisons, list) or len(comparisons) > 20:
        raise Invalid('bounded comparisons required')
    if window_start is not None:
        day(window_start, 'window_start')
    selected = selection(dataset, view, stage, window_start)
    episodes = selected['episodes']
    curve = cumulative_incidence(_subject(episodes))
    results = [comparison(episodes, spec, target) for spec in comparisons]
    ids = [row['id'] for row in results]
    if len(ids) != len(set(ids)):
        raise Invalid('duplicate comparison id')
    return {'stage': stage, 'target_day': target,
            'definition': {'stage': stage,
                           'window_start': selected['window']['start'],
                           'window_end': view['as_of'],
                           'as_of': view['as_of'], 'known_at': view['known_at'],
                           'target_day': target,
                           'time_origin': 'admission',
                           'outcomes': list(OUTCOMES),
                           'selection': 'one_first_visible_inpatient_episode_per_person',
                           'censoring': 'open_at_business_cutoff',
                           'ties': 'events_before_censoring',
                           'formal_estimand': 'aalen_johansen_cumulative_incidence',
                           'comparisons': comparisons},
            'selected': selected,
            'episodes': episodes,
            'curve': curve,
            'target': target_at(curve, target),
            'aligned_naive': aligned_naive(_subject(episodes), curve),
            'naive': naive_estimates(_subject(episodes), target),
            'comparisons': results,
            'notes': ['descriptive_not_clinical_recommendation',
                      'competing_discharges_are_not_censoring',
                      'naive_estimates_shown_only_to_explain_difference',
                      'target_beyond_supported_followup_is_not_estimable']}


def verify_study(dataset, view, stage, result):
    if result.get('stage') != stage:
        raise Invalid('competing-risk stage mismatch')
    verify(result['curve'])
    selected = result['selected']
    episodes = selected['episodes']
    if len({row['person'] for row in episodes}) != len(episodes):
        raise Invalid('competing-risk cohort repeats a person')
    by_person = defaultdict(list)
    for key, row in view['visits'].items():
        if dataset.stage(key) == stage:
            by_person[row['person']].append(key)
    covered = set()
    for row in episodes:
        key = row['visit']
        source = view['visits'].get(key)
        if source is None or key not in by_person[row['person']]:
            raise Invalid('competing-risk episode missing from visible source')
        if row['source_digest'] != digest(source):
            raise Invalid('competing-risk episode source digest mismatch')
        if row['admitted'] != source['admitted'] or row['discharged'] != source['discharged']:
            raise Invalid('competing-risk episode dates differ from source')
        expected_duration = days(source['admitted'], source['discharged'] or view['as_of'])
        if row['duration'] != expected_duration:
            raise Invalid('competing-risk duration differs from source')
        if row['outcome'] != (source['outcome'] if source['discharged'] else None):
            raise Invalid('competing-risk outcome differs from source')
        if covered & {key}:
            raise Invalid('competing-risk episode selected twice')
        covered.add(key)
    excluded = [row['visit'] for row in selected['excluded_records']]
    invisible = [row['visit'] for row in selected['invisible_records']]
    if len(excluded) != len(set(excluded)) or set(covered) & set(excluded):
        raise Invalid('competing-risk exclusion listing overlaps selected episodes')
    if set(invisible) - set(excluded):
        raise Invalid('invisible competing-risk records are not listed as excluded')
    visible_stage = {key for key, row in view['visits'].items() if dataset.stage(key) == stage}
    if covered | set(excluded) != visible_stage | set(invisible):
        raise Invalid('competing-risk selection does not partition stage records')
    if result['target']['day'] != result['definition']['target_day']:
        raise Invalid('competing-risk target day mismatch')
    expected = target_at(result['curve'], result['target']['day'])
    if result['target'] != expected:
        raise Invalid('competing-risk target estimate mismatch')
    for row in result['comparisons']:
        spec = next((item for item in result['definition']['comparisons'] if item['id'] == row['id']), None)
        if spec is None or row['outcome'] != spec['outcome']:
            raise Invalid('competing-risk comparison missing from definition')
        for side in ('left', 'right'):
            if not row[side]['estimable'] and row[side]['cumulative'] is not None:
                raise Invalid('unsupported comparison cohort has an estimate')
        bootstrap = row['bootstrap']
        if bootstrap['supported_replicates'] + bootstrap['unsupported_replicates'] != bootstrap['replicates']:
            raise Invalid('bootstrap did not retain unsupported replicates')
        if bootstrap['supported_replicates'] and bootstrap['interval'] is None:
            raise Invalid('supported bootstrap replicates lack an interval')
    return {'valid': True, 'episodes': len(episodes), 'comparisons': len(result['comparisons'])}

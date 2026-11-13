from encounters.values import digest, days, Invalid
from pathways.observations import referral_waits
from pathways.chains import build


def waiting_lineage(dataset, report, referral_id):
    options = report['options']
    view = dataset.slice(options['as_of'], options['known_at'])
    if referral_id not in view['referrals']:
        raise KeyError(referral_id)
    referral = view['referrals'][referral_id]
    source = view['visits'][referral['source']]
    target = view['visits'].get(referral['target'])
    observation = next((row for row in report['wait_observations']['rows'] if row['id'] == referral_id), None)
    excluded = next((row for row in report['wait_observations']['excluded'] if row['id'] == referral_id), None)
    records = [{'table': 'referrals', 'id': referral_id, 'digest': digest(referral)},
               {'table': 'visits', 'id': source['id'], 'digest': digest(source)},
               {'table': 'people', 'id': referral['person'], 'digest': digest(dataset.people[referral['person']])},
               {'table': 'facilities', 'id': source['facility'], 'digest': digest(dataset.facilities[source['facility']])}]
    if target:
        records.append({'table': 'visits', 'id': target['id'], 'digest': digest(target)})
    return {'referral': referral_id, 'records': records, 'observation': observation,
            'excluded': excluded, 'source_cutoffs': options,
            'formula': 'target_admission_minus_request' if target else 'as_of_minus_request',
            'unit_selection': report['waiting']['unit'],
            'excluded_as_repeated_person': referral_id in report['waiting']['excluded_repeat_referrals']}


def chain_lineage(dataset, report, root):
    chain = next((row for row in report['chains']['chains'] if row['id'] == root), None)
    if chain is None:
        raise KeyError(root)
    view = dataset.slice(report['options']['as_of'], report['options']['known_at'])
    records = [{'table': 'visits', 'id': key, 'digest': digest(view['visits'][key])} for key in chain['visits']]
    for edge in chain['edges']:
        records.append({'table': 'referrals', 'id': edge['referral'],
                        'digest': digest(view['referrals'][edge['referral']])})
    return {'root': root, 'chain_digest': chain['chain_digest'], 'records': records,
            'care_day_terms': [{'visit': key, 'days': days(view['visits'][key]['admitted'],
                               view['visits'][key]['discharged'] or view['as_of'])} for key in chain['visits']],
            'gap_day_terms': [{'referral': edge['referral'], 'days': edge['handoff_gap_days']} for edge in chain['edges']],
            'elapsed_days': chain['elapsed_days']}


def evidence_index(dataset, report):
    waits = []
    for row in report['wait_observations']['rows']:
        evidence = waiting_lineage(dataset, report, row['id'])
        waits.append({'referral': row['id'], 'lineage_digest': digest(evidence),
                      'included_in_selected_unit': not evidence['excluded_as_repeated_person']})
    chains = []
    for chain in report['chains']['chains']:
        evidence = chain_lineage(dataset, report, chain['id'])
        chains.append({'root': chain['id'], 'lineage_digest': digest(evidence)})
    body = {'input_digest': dataset.fingerprint, 'report_digest': report['report_digest'],
            'waiting_evidence': waits, 'chain_evidence': chains}
    return {**body, 'evidence_digest': digest(body)}


def verify_index(dataset, report, index):
    expected = evidence_index(dataset, report)
    if index != expected:
        raise Invalid('provenance index does not match source replay')
    return {'valid': True, 'waiting': len(index['waiting_evidence']), 'chains': len(index['chain_evidence'])}


def competing_lineage(dataset, report, visit_id):
    study = report['competing_risks']
    if visit_id not in dataset.visits:
        raise KeyError(visit_id)
    view = dataset.slice(report['options']['as_of'], report['options']['known_at'])
    episode = next((row for row in study['episodes'] if row['visit'] == visit_id), None)
    selected = episode is not None
    excluded = None if selected else next(
        (row for row in study['selected']['excluded_records'] if row['visit'] == visit_id), None)
    invisible = next((row for row in study['selected']['invisible_records'] if row['visit'] == visit_id), None)
    outside_stage = dataset.stage(visit_id) != study['stage']
    records = []
    if visit_id in view['visits']:
        source = view['visits'][visit_id]
        records.append({'table': 'visits', 'id': visit_id, 'digest': digest(source)})
        records.append({'table': 'people', 'id': source['person'],
                        'digest': digest(dataset.people[source['person']])})
        records.append({'table': 'facilities', 'id': source['facility'],
                        'digest': digest(dataset.facilities[source['facility']])})
    if episode is not None and visit_id not in view['visits']:
        raise Invalid('selected competing-risk episode is not visible')
    term = None
    if episode is not None:
        term = {'visit': visit_id, 'admitted': episode['admitted'],
                'end': episode['discharged'] or view['as_of'],
                'duration': episode['duration'],
                'outcome': episode['outcome'],
                'censoring_reason': episode['censoring_reason'],
                'formula': 'discharge_or_as_of_minus_admission'}
        curve = study['curve']
        terms = []
        for row in curve['rows']:
            for outcome in curve['outcomes']:
                count = row['outcome_events'][outcome]
                if count:
                    terms.append({'time': row['time'], 'outcome': outcome, 'events': count,
                                  'at_risk': row['at_risk'],
                                  'increment': {'numerator': count, 'denominator': row['at_risk']}})
        term['cumulative_incidence_terms_for_episode'] = [
            row for row in terms if row['time'] == episode['duration']
            and row['outcome'] == episode['outcome']]
        term['risk_set_terms_at_episode_time'] = [row for row in terms
                                                  if row['time'] == episode['duration']]
    return {'visit': visit_id, 'selected': selected, 'episode': episode,
            'excluded': excluded, 'invisible': invisible, 'records': records,
            'term': term, 'stage': study['stage'], 'outside_selected_stage': outside_stage,
            'definition': study['definition'], 'counts_to_curve': selected,
            'note': 'excluded_or_invisible_records_do_not_enter_any_risk_set'}

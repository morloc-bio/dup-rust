// Rust implementations of biodup.

#[derive(Clone, Debug, Default)]
pub struct DupKeys {
    pub keyed: Vec<bool>,
    pub words: Vec<u64>,
}

#[allow(non_snake_case)]
#[derive(Clone, Debug, Default)]
pub struct DupRuns {
    pub first: Vec<u64>,
    pub firstCount: i64,
    pub last: Vec<u64>,
    pub lastCount: i64,
    pub interior: Vec<i64>,
}

fn dup_word(kind: u64, plus: bool, seq: i64, pos: i64) -> u64 {
    if !(0..1 << 29).contains(&seq) {
        panic!("duplicate keys: genome sequence index {} is outside 0 .. 2^29", seq);
    }
    if !(0..1 << 32).contains(&pos) {
        panic!("duplicate keys: genome position {} is outside 0 .. 2^32", pos);
    }
    kind << 62 | (plus as u64) << 61 | (seq as u64) << 32 | pos as u64
}

fn dup_check_pos(pos: i64) -> u64 {
    if !(0..1 << 32).contains(&pos) {
        panic!("duplicate keys: genome position {} is outside 0 .. 2^32", pos);
    }
    pos as u64
}

// A fragment's key, from its best placement and its mates' lengths.
fn dup_key(refs: &BioRefModels, f: &BioFragHits, l1: i64, l2: i64) -> Option<[u64; 2]> {
    let mut best: Option<&BioTxHit> = None;
    for h in &f.hits {
        if best.map_or(true, |b| h.weight > b.weight) {
            best = Some(h);
        }
    }
    let h = best?;
    let m = refs.maps.get(h.tx as usize)?;
    if m.dst < 0 {
        return None;
    }
    if h.mates == 0 && h.fragLen > 0 {
        let (a, sa) = bio_project_point(m, h.pos)?;
        let (b, _) = bio_project_point(m, h.pos + h.fragLen - 1)?;
        let read1_plus = h.fwd == sa;
        return Some([dup_word(0, read1_plus, m.dst, a.min(b)), dup_check_pos(a.max(b))]);
    }
    let (mate, len) = match h.mates {
        2 => (2, l2),
        1 => (1, l1),
        _ => (0, l1),
    };
    let five = if h.fwd { h.pos } else { h.pos + len - 1 };
    let (g, s) = bio_project_point(m, five)?;
    Some([dup_word(1 + mate, h.fwd == s, m.dst, g), 0])
}

fn dup_keys(refs: &BioRefModels, lens: &[(i64, i64)], hits: &[BioFragHits]) -> DupKeys {
    if lens.len() != hits.len() {
        panic!("dupKeys: {} reads but {} placements", lens.len(), hits.len());
    }
    let mut out = DupKeys { keyed: Vec::with_capacity(hits.len()), words: Vec::with_capacity(2 * hits.len()) };
    for (f, &(a, b)) in hits.iter().zip(lens) {
        match dup_key(refs, f, a, b) {
            Some(k) => {
                out.keyed.push(true);
                out.words.extend_from_slice(&k);
            }
            None => out.keyed.push(false),
        }
    }
    out
}

pub fn dup_keys_pairs(_anchor: i64, refs: &BioRefModels, reads: &Vec<(BioFastqEntry, BioFastqEntry)>, hits: &Vec<BioFragHits>) -> DupKeys {
    let lens: Vec<(i64, i64)> = reads.iter().map(|(a, b)| (a.bases.bases.len() as i64, b.bases.bases.len() as i64)).collect();
    dup_keys(refs, &lens, hits)
}

pub fn dup_keys_reads(_anchor: i64, refs: &BioRefModels, reads: &Vec<BioFastqEntry>, hits: &Vec<BioFragHits>) -> DupKeys {
    let lens: Vec<(i64, i64)> = reads.iter().map(|a| (a.bases.bases.len() as i64, 0)).collect();
    dup_keys(refs, &lens, hits)
}

fn dup_add(h: &mut Vec<i64>, n: i64, k: i64) {
    let n = n as usize;
    if h.len() <= n {
        h.resize(n + 1, 0);
    }
    h[n] += k;
}

// The runs of a summary in order: (key, length), the interior as a histogram.
fn dup_ends(r: &DupRuns) -> (Option<(&[u64], i64)>, Option<(&[u64], i64)>) {
    if r.first.is_empty() {
        return (None, None);
    }
    let head = Some((r.first.as_slice(), r.firstCount));
    if r.last.is_empty() {
        (head, None)
    } else {
        (head, Some((r.last.as_slice(), r.lastCount)))
    }
}

pub fn dup_runs(xs: &Vec<u64>) -> DupRuns {
    if xs.len() % 2 != 0 {
        panic!("dupRuns: {} words are not whole keys", xs.len());
    }
    let mut runs: Vec<(&[u64], i64)> = Vec::new();
    for k in xs.chunks_exact(2) {
        match runs.last_mut() {
            Some((p, n)) if *p == k => *n += 1,
            Some((p, _)) if *p > k => panic!("dupRuns: keys are not ascending"),
            _ => runs.push((k, 1)),
        }
    }
    let mut out = DupRuns::default();
    if let Some(&(k, n)) = runs.first() {
        out.first = k.to_vec();
        out.firstCount = n;
    }
    if runs.len() > 1 {
        let &(k, n) = runs.last().unwrap();
        out.last = k.to_vec();
        out.lastCount = n;
        for &(_, n) in &runs[1..runs.len() - 1] {
            dup_add(&mut out.interior, n, 1);
        }
    }
    out
}

pub fn dup_merge_runs(a: &DupRuns, b: &DupRuns) -> DupRuns {
    let (a_head, a_last) = dup_ends(a);
    let (b_head, b_last) = dup_ends(b);
    let Some(a_head) = a_head else { return b.clone() };
    let Some(b_head) = b_head else { return a.clone() };
    let mut interior = a.interior.clone();
    for (n, &k) in b.interior.iter().enumerate() {
        if k > 0 {
            dup_add(&mut interior, n as i64, k);
        }
    }
    // The runs in order: a_head, [a interior], a_last?, b_head, [b interior], b_last?
    let a_tail = a_last.unwrap_or(a_head);
    if a_tail.0 > b_head.0 {
        panic!("mergeDupRuns: batches are not in ascending order");
    }
    let joined = a_tail.0 == b_head.0;
    let mut first = a_head;
    let mut middle: Vec<i64> = Vec::new();
    let mut last: Option<(&[u64], i64)>;
    if joined {
        let j = (a_tail.0, a_tail.1 + b_head.1);
        match (a_last, b_last) {
            (None, None) => {
                first = j;
                last = None;
            }
            (None, Some(bl)) => {
                first = j;
                last = Some(bl);
            }
            (Some(_), None) => {
                last = Some(j);
            }
            (Some(_), Some(bl)) => {
                middle.push(j.1);
                last = Some(bl);
            }
        }
    } else {
        if let Some(al) = a_last {
            middle.push(al.1);
        }
        match b_last {
            None => last = Some(b_head),
            Some(bl) => {
                middle.push(b_head.1);
                last = Some(bl);
            }
        }
    }
    for n in middle {
        dup_add(&mut interior, n, 1);
    }
    let mut out = DupRuns { first: first.0.to_vec(), firstCount: first.1, last: Vec::new(), lastCount: 0, interior };
    if let Some(l) = last.take() {
        out.last = l.0.to_vec();
        out.lastCount = l.1;
    }
    out
}

pub fn dup_histogram(r: &DupRuns) -> Vec<i64> {
    let mut h = r.interior.clone();
    let (head, last) = dup_ends(r);
    for (_, n) in head.into_iter().chain(last) {
        dup_add(&mut h, n, 1);
    }
    while h.len() > 1 && *h.last().unwrap() == 0 {
        h.pop();
    }
    if h.is_empty() {
        h.push(0);
    }
    h
}

fn dup_totals(h: &[i64]) -> (i64, i64) {
    let distinct: i64 = h.iter().skip(1).sum();
    let total: i64 = h.iter().enumerate().map(|(j, &k)| j as i64 * k).sum();
    (distinct, total)
}

pub fn dup_fraction(h: &Vec<i64>) -> f64 {
    let (distinct, total) = dup_totals(h);
    if total == 0 { 0.0 } else { 1.0 - distinct as f64 / total as f64 }
}

// The probability that none of the j copies of a key is among n fragments
// drawn without replacement from `total`.
fn dup_missed(total: i64, n: i64, j: i64) -> f64 {
    if n + j > total {
        return 0.0;
    }
    let mut p = 1.0f64;
    for i in 0..j {
        p *= (total - n - i) as f64 / (total - i) as f64;
        if p < 1e-300 {
            return 0.0;
        }
    }
    p
}

pub fn dup_complexity_curve(step: i64, h: &Vec<i64>) -> Vec<(i64, f64)> {
    if step <= 0 {
        panic!("complexityCurve: step must be positive, not {}", step);
    }
    let (_, total) = dup_totals(h);
    let mut out = Vec::new();
    let mut n = step.min(total);
    while total > 0 {
        let e: f64 = h.iter().enumerate().skip(1).filter(|(_, &k)| k > 0).map(|(j, &k)| k as f64 * (1.0 - dup_missed(total, n, j as i64))).sum();
        out.push((n, e));
        if n == total {
            break;
        }
        n = (n + step).min(total);
    }
    out
}

pub fn dup_rate_text(h: &Vec<i64>) -> String {
    let mut s = String::from("Occurrence\tUniqReadNumber\n");
    for (j, &k) in h.iter().enumerate().skip(1) {
        if k > 0 {
            s.push_str(&format!("{}\t{}\n", j, k));
        }
    }
    s
}

pub fn dup_complexity_curve_text(step: i64, h: &Vec<i64>) -> String {
    let mut s = String::from("TOTAL_READS\tEXPECTED_DISTINCT\n");
    for (n, e) in dup_complexity_curve(step, h) {
        s.push_str(&format!("{}\t{:.1}\n", n, e));
    }
    s
}

pub fn dup_rust_anchor(x: i64) -> i64 {
    x
}

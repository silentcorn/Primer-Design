from primer import temp_diff, is_length, runs_repeats

def test_temp_diff():
    assert temp_diff(54, 50) == 0

def test_is_length():
    assert is_length('123456789') == 2

def test_runs_repeats():
    seq = 'GCGCGCGCGCGCGCGCGC'
    assert runs_repeats(seq) == 1
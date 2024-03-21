import pytest
from program.unknown_bases_stats import UnknownBasesStats

class GzipFileObject:
    _LINES = []
    def __init__(self, _1, _2) -> None:
        pass
    
    def __enter__(self):
        return GzipFileObject._LINES
    
    def __exit__(self, _1, _2, _3):
        return

def test_run_continuing_across_lines_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar\n", f"{'N'*5}\n", f"{'N'*5}\n"]
    # Act
    sut = UnknownBasesStats('foo')
    result = sut._count_unknown_bases()
    # Assert
    assert '1' in result, f"Expected sequence with name '1' but not found"
    assert len(result) == 1, f"Expected 1 sequence but got {len(result)}"
    assert len(result['1']) == 1, f"Expected 1 run of Ns but got {len(result['1'])}"
    assert result['1'][0][0] == 0, f"Expected 1st run to start at 0 but it starts at {result['1'][0][0]}"
    assert result['1'][0][1] == 10, f"Expected 1st run to end at 10 but it ends at {result['1'][0][1]}"
    
def test_run_starting_at_0_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar\n", f"{'N'*3}{'A'*5}\n"]
    # Act
    sut = UnknownBasesStats('foo')
    result = sut._count_unknown_bases()
    # Assert
    assert '1' in result, f"Expected sequence with name '1' but not found"
    assert len(result) == 1, f"Expected 1 sequence but got {len(result)}"
    assert len(result['1']) == 1, f"Expected 1 run of Ns but got {len(result['1'])}"
    assert result['1'][0][0] == 0, f"Expected 1st run to start at 0 but it starts at {result['1'][0][0]}"
    assert result['1'][0][1] == 3, f"Expected 1st run lenght to be 3 but it is {result['1'][0][1]}"
    
def test_run_starting_in_the_middle_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar\n", f"{'A'*3}{'N'*5}{'A'*3}\n"]
    # Act
    sut = UnknownBasesStats('foo')
    result = sut._count_unknown_bases()
    # Assert
    assert '1' in result, f"Expected sequence with name '1' but not found"
    assert len(result) == 1, f"Expected 1 sequence but got {len(result)}"
    assert len(result['1']) == 1, f"Expected 1 run of Ns but got {len(result['1'])}"
    assert result['1'][0][0] == 3, f"Expected 1st run to start at 3 but it starts at {result['1'][0][0]}"
    assert result['1'][0][1] == 5, f"Expected 1st run lenght to be 5 but it is {result['1'][0][1]}"

def test_run_ending_with_line_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar", f"{'A'*3}{'N'*5}\n"]
    # Act
    sut = UnknownBasesStats('foo')
    result = sut._count_unknown_bases()
    # Assert
    assert '1' in result, f"Expected sequence with name '1' but not found"
    assert len(result) == 1, f"Expected 1 sequence but got {len(result)}"
    assert len(result['1']) == 1, f"Expected 1 runs of Ns but got {len(result['1'])}"
    assert result['1'][0][0] == 3, f"Expected 1st run to start at 3 but it starts at {result['1'][0][0]}"
    assert result['1'][0][1] == 5, f"Expected 1st run lenght to be 5 but it is {result['1'][0][1]}"
    
def test_run_not_continuing_in_next_line_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar", f"{'A'*3}{'N'*5}\n", f"{'A'*3}\n"]
    # Act
    sut = UnknownBasesStats('foo')
    result = sut._count_unknown_bases()
    # Assert
    assert '1' in result, f"Expected sequence with name '1' but not found"
    assert len(result) == 1, f"Expected 1 sequence but got {len(result)}"
    assert len(result['1']) == 1, f"Expected 1 runs of Ns but got {len(result['1'])}"
    assert result['1'][0][0] == 3, f"Expected 1st run to start at 3 but it starts at {result['1'][0][0]}"
    assert result['1'][0][1] == 5, f"Expected 1st run lenght to be 5 but it is {result['1'][0][1]}"
    
def test_run_not_continuing_in_next_line_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar", f"{'A'*3}{'N'*5}\n", ">2 foo bar", f"{'A'*3}{'N'*5}\n"]
    # Act
    sut = UnknownBasesStats('foo')
    result = sut._count_unknown_bases()
    # Assert
    assert '1' in result, f"Expected sequence with name '1' but not found"
    assert '2' in result, f"Expected sequence with name '2' but not found"
    assert len(result) == 2, f"Expected 2 sequence but got {len(result)}"
    assert len(result['1']) == 1, f"Expected 1 run in the 1st sequence but got {len(result['1'])}"
    assert result['1'][0][0] == 3, f"Expected 1st run of 1st sequence to start at 3 but it starts at {result['1'][0][0]}"
    assert result['1'][0][1] == 5, f"Expected 1st run of 1st sequence lenght to be 5 but it is {result['1'][0][1]}"
    assert len(result['2']) == 1, f"Expected 1 run in the 2nd sequence but got {len(result['2'])}"
    assert result['2'][0][0] == 3, f"Expected 1st run of 2nd sequence to start at 3 but it starts at {result['2'][0][0]}"
    assert result['2'][0][1] == 5, f"Expected 1st run of 2nd sequence lenght to be 5 but it is {result['2'][0][1]}"
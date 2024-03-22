from pathlib import Path
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
    sut = UnknownBasesStats(Path('foo'))
    result = sut._count_unknown_bases()
    # Assert
    assert len([x for x in result if x.name is '1']) == 1, f"Expected sequence with name '1' but not found"
    assert len(result[0].runs) == 1, f"Expected 1 run of Ns but got {len(result[0].runs)}"
    assert result[0].runs[0].start == 0, f"Expected 1st run to start at 0 but it starts at {result[0].runs[0].start}"
    assert result[0].runs[0].lenght == 10, f"Expected 1st run to end at 10 but it ends at {result[0].runs[0].lenght}"
    
def test_run_starting_at_0_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar\n", f"{'N'*3}{'A'*5}\n"]
    # Act
    sut = UnknownBasesStats(Path('foo'))
    result = sut._count_unknown_bases()
    # Assert
    assert len([x for x in result if x.name is '1']) == 1, f"Expected sequence with name '1' but not found"
    assert len(result[0].runs) == 1, f"Expected 1 run of Ns but got {len(result[0].runs)}"
    assert result[0].runs[0].start == 0, f"Expected 1st run to start at 0 but it starts at {result[0].runs[0].start}"
    assert result[0].runs[0].lenght == 3, f"Expected 1st run lenght to be 3 but it is {result[0].runs[0].lenght}"
    
def test_run_starting_in_the_middle_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar\n", f"{'A'*3}{'N'*5}{'A'*3}\n"]
    # Act
    sut = UnknownBasesStats(Path('foo'))
    result = sut._count_unknown_bases()
    # Assert
    assert len([x for x in result if x.name is '1']) == 1, f"Expected sequence with name '1' but not found"
    assert len(result[0].runs) == 1, f"Expected 1 run of Ns but got {len(result[0].runs)}"
    assert result[0].runs[0].start == 3, f"Expected 1st run to start at 3 but it starts at {result[0].runs[0].start}"
    assert result[0].runs[0].lenght == 5, f"Expected 1st run lenght to be 5 but it is {result[0].runs[0].lenght}"

def test_run_ending_with_line_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar", f"{'A'*3}{'N'*5}\n"]
    # Act
    sut = UnknownBasesStats(Path('foo'))
    result = sut._count_unknown_bases()
    # Assert
    assert len([x for x in result if x.name is '1']) == 1, f"Expected sequence with name '1' but not found"
    assert len(result[0].runs) == 1, f"Expected 1 runs of Ns but got {len(result[0].runs)}"
    assert result[0].runs[0].start == 3, f"Expected 1st run to start at 3 but it starts at {result[0].runs[0].start}"
    assert result[0].runs[0].lenght == 5, f"Expected 1st run lenght to be 5 but it is {result[0].runs[0].lenght}"
    
def test_run_not_continuing_in_next_line_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar", f"{'A'*3}{'N'*5}\n", f"{'A'*3}\n"]
    # Act
    sut = UnknownBasesStats(Path('foo'))
    result = sut._count_unknown_bases()
    # Assert
    assert len([x for x in result if x.name is '1']) == 1, f"Expected sequence with name '1' but not found"
    assert len(result[0].runs) == 1, f"Expected 1 runs of Ns but got {len(result[0].runs)}"
    assert result[0].runs[0].start == 3, f"Expected 1st run to start at 3 but it starts at {result[0].runs[0].start}"
    assert result[0].runs[0].lenght == 5, f"Expected 1st run lenght to be 5 but it is {result[0].runs[0].lenght}"
    
def test_run_not_continuing_in_next_line_is_processed_correctly():
    # Arrange
    pytest.MonkeyPatch().setattr('gzip.open', GzipFileObject)
    GzipFileObject._LINES = [">1 foo bar", f"{'A'*3}{'N'*5}\n", ">2 foo bar", f"{'A'*3}{'N'*5}\n"]
    # Act
    sut = UnknownBasesStats(Path('foo'))
    result = sut._count_unknown_bases()
    # Assert
    assert len([x for x in result if x.name is '1']) == 1, f"Expected sequence with name '1' but not found"
    assert len([x for x in result if x.name is '2']) == 1, f"Expected sequence with name '2' but not found"
    assert len(result[0].runs) == 1, f"Expected 1 run in the 1st sequence but got {len(result[0].runs)}"
    assert result[0].runs[0].start == 3, f"Expected 1st run of 1st sequence to start at 3 but it starts at {result[0].runs[0].start}"
    assert result[0].runs[0].lenght == 5, f"Expected 1st run of 1st sequence lenght to be 5 but it is {result[0].runs[0].lenght}"
    assert len(result[1].runs) == 1, f"Expected 1 run in the 2nd sequence but got {len(result[1].runs)}"
    assert result[1].runs[0].start == 3, f"Expected 1st run of 2nd sequence to start at 3 but it starts at {result[1].runs[0].start}"
    assert result[1].runs[0].lenght == 5, f"Expected 1st run of 2nd sequence lenght to be 5 but it is {result[1].runs[0].lenght}"
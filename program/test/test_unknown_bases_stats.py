import pytest
from program.reference.fasta_file import FastaFile


class MockFile:
    def __init__(self, lines) -> None:
        self.lines = lines

    def __enter__(self):
        return self.lines

    def __exit__(self, _1, _2, _3):
        return


class MockPath:
    def __init__(self, lines) -> None:
        self.lines = lines
        self.name = "foo.fa.gz"
        self.stem = "foo.fa"

    def open(self, *_1):
        return MockFile(self.lines)

    def exists(self):
        return True


class MockGenome:
    def __init__(self, fasta_lines, dictionary_lines) -> None:
        self.final_name = MockPath(fasta_lines)
        self.dict = MockPath(dictionary_lines)


def gzip_open(path, mode):
    return MockFile(path.lines)


def test_run_continuing_across_lines_is_processed_correctly():
    # Arrange
    fa_lines = [
        ">1 irrelevant\n",
        f"{'N'*5}\n",
        f"{'N'*5}\n",
    ]
    dict_line = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:10 M5:dummy UR:file://c:/foo.fa.gz",
    ]
    pytest.MonkeyPatch().setattr("gzip.open", gzip_open)
    sut = FastaFile(MockGenome(fa_lines, dict_line))

    # Act
    result = sut.count_letters()

    # Assert
    assert len([x for x in result if x.name == "1"]) == 1
    assert len(result[0].runs) == 1
    assert result[0].runs[0].start == 0
    assert result[0].runs[0].length == 10


def test_run_starting_at_0_is_processed_correctly():
    # Arrange
    fa_lines = [">1 irrelevant\n", f"{'N'*3}{'A'*5}\n"]
    dict_line = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:8 M5:dummy UR:file://c:/foo.fa.gz",
    ]
    pytest.MonkeyPatch().setattr("gzip.open", gzip_open)
    sut = FastaFile(MockGenome(fa_lines, dict_line))

    # Act
    result = sut.count_letters()

    # Assert
    assert len([x for x in result if x.name == "1"]) == 1
    assert len(result[0].runs) == 1
    assert result[0].runs[0].start == 0
    assert result[0].runs[0].length == 3


def test_run_starting_in_the_middle_is_processed_correctly():
    # Arrange
    fa_lines = [
        ">1 irrelevant\n",
        f"{'A'*3}{'N'*5}{'A'*3}\n",
    ]
    dict_line = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:11 M5:dummy UR:file://c:/foo.fa.gz",
    ]
    pytest.MonkeyPatch().setattr("gzip.open", gzip_open)
    sut = FastaFile(MockGenome(fa_lines, dict_line))

    # Act
    result = sut.count_letters()

    # Assert
    assert len([x for x in result if x.name == "1"]) == 1
    assert len(result[0].runs) == 1
    assert result[0].runs[0].start == 3
    assert result[0].runs[0].length == 5


def test_run_ending_with_line_is_processed_correctly():
    # Arrange
    fa_lines = [">1 irrelevant", f"{'A'*3}{'N'*5}\n"]
    dict_line = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:8 M5:dummy UR:file://c:/foo.fa.gz",
    ]
    pytest.MonkeyPatch().setattr("gzip.open", gzip_open)
    sut = FastaFile(MockGenome(fa_lines, dict_line))

    # Act
    result = sut.count_letters()

    # Assert
    assert len([x for x in result if x.name == "1"]) == 1
    assert len(result[0].runs) == 1
    assert result[0].runs[0].start == 3
    assert result[0].runs[0].length == 5


def test_run_not_continuing_in_next_line_is_processed_correctly():
    # Arrange
    fa_lines = [
        ">1 irrelevant",
        f"{'A'*3}{'N'*5}\n",
        f"{'A'*3}\n",
    ]
    dict_line = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:11 M5:dummy UR:file://c:/foo.fa.gz",
    ]
    pytest.MonkeyPatch().setattr("gzip.open", gzip_open)
    sut = FastaFile(MockGenome(fa_lines, dict_line))

    # Act
    result = sut.count_letters()

    # Assert
    assert len([x for x in result if x.name == "1"]) == 1
    assert len(result[0].runs) == 1
    assert result[0].runs[0].start == 3
    assert result[0].runs[0].length == 5


def test_run_not_continuing_in_next_line_is_processed_correctly():
    # Arrange
    fa_lines = [
        ">1 irrelevant",
        f"{'A'*3}{'N'*5}\n",
        ">2 irrelevant",
        f"{'A'*3}{'N'*5}\n",
    ]
    dict_line = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:8 M5:dummy UR:file://c:/foo.fa.gz",
        "@SQ SN:2 LN:8 M5:dummy UR:file://c:/foo.fa.gz",
    ]
    pytest.MonkeyPatch().setattr("gzip.open", gzip_open)
    sut = FastaFile(MockGenome(fa_lines, dict_line))

    # Act
    result = sut.count_letters()

    # Assert
    assert len([x for x in result if x.name == "1"]) == 1
    assert len([x for x in result if x.name == "2"]) == 1
    assert len(result[0].runs) == 1
    assert result[0].runs[0].start == 3
    assert result[0].runs[0].length == 5
    assert len(result[1].runs) == 1
    assert result[1].runs[0].start == 3
    assert result[1].runs[0].length == 5

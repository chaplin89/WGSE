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
    def __init__(self, lines, exists=True) -> None:
        self.lines = lines
        self.name = "foo.fa.gz"
        self.stem = "foo.fa"
        self._exists = exists

    def open(self, *_1):
        return MockFile(self.lines)

    def exists(self):
        return self._exists


class MockGenome:
    def __init__(self, fasta, dict) -> None:
        self.final_name = fasta
        self.dict = dict


def gzip_open(path, mode):
    return MockFile(path.lines)


def test_no_dictionary():
    fa_lines = [
        ">1 irrelevant\n",
        f"{'N'*5}\n",
        f"{'N'*5}\n",
        ">2 irrelevant\n",
    ]
    dict_lines = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:10 M5:dummy UR:file://c:/foo.fa.gz",
    ]
   
    with pytest.raises(RuntimeError) as e:
        genome = MockGenome(MockPath(fa_lines), MockPath(None,False))
        sut = FastaFile(genome)
    assert "Unable to find dictionary" in str(e.value)

def test_sequence_not_in_dictionary():
    fa_lines = [
        ">1 irrelevant\n",
        f"{'N'*5}\n",
        f"{'N'*5}\n",
        ">2 irrelevant\n",
    ]
    dict_lines = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:10 M5:dummy UR:file://c:/foo.fa.gz",
    ]
    pytest.MonkeyPatch().setattr("gzip.open", gzip_open)
   
    with pytest.raises(ValueError) as e:
        genome = MockGenome(MockPath(fa_lines), MockPath(dict_lines))
        FastaFile(genome).split_into_sequences()
    assert "not present in dictionary" in str(e.value)
    
def test_fastq():
    fa_lines = [
        ">1 irrelevant\n",
        f"{'N'*5}\n",
        f"{'N'*5}\n",
        "+2 irrelevant\n",
    ]
    dict_lines = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:10 M5:dummy UR:file://c:/foo.fa.gz",
    ]
    pytest.MonkeyPatch().setattr("gzip.open", gzip_open)
   
    with pytest.raises(RuntimeError) as e:
        genome = MockGenome(MockPath(fa_lines), MockPath(dict_lines))
        FastaFile(genome).split_into_sequences()
    assert "Expected a FASTA" in str(e.value)
    
def test_duplicate_sequence():
    fa_lines = [
        ">1 irrelevant\n",
        f"{'N'*5}\n",
        f"{'N'*5}\n",
        ">1 irrelevant\n",
    ]
    dict_lines = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:10 M5:dummy UR:file://c:/foo.fa.gz",
    ]
    pytest.MonkeyPatch().setattr("gzip.open", gzip_open)
   
    with pytest.raises(RuntimeError) as e:
        genome = MockGenome(MockPath(fa_lines), MockPath(dict_lines))
        FastaFile(genome).split_into_sequences()
    assert "duplicated sequence" in str(e.value)
    
    
def test_only_comments():
    fa_lines = [
        "#Hello\n",
        f"#Foo\n"
    ]
    dict_lines = [
        "@HD VN:1.0 SO:unsorted",
        "@SQ SN:1 LN:10 M5:dummy UR:file://c:/foo.fa.gz",
    ]
    pytest.MonkeyPatch().setattr("gzip.open", gzip_open)
   
    with pytest.raises(RuntimeError) as e:
        genome = MockGenome(MockPath(fa_lines), MockPath(dict_lines))
        FastaFile(genome).split_into_sequences()
    assert "no sequences" in str(e.value)
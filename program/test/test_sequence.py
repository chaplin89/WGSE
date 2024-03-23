import pytest
from program.unknown_bases_stats import Sequence


def test_sequence_with_invalid_length():
    with pytest.raises(IndexError) as e:
        sut = Sequence("Foo", 0)
    assert "greater than" in str(e.value)

def test_close_run_less_equal_than_open():
    with pytest.raises(RuntimeError) as e:
        sut = Sequence("Foo", 10)
        sut.open_run(0)
        sut.close_run(0)
    assert "with position" in str(e.value)

def test_close_run_greater_than_open():
    # Arrange
    sut = Sequence("Foo", 10)
    # Act    
    sut.open_run(0)
    sut.close_run(1)
    # Assert
    assert len(sut.runs) == 1
    assert sut.runs[0].start == 0
    assert sut.runs[0].length == 1
        
def test_double_open_run():
    with pytest.raises(RuntimeError) as e:
        sut = Sequence("Foo", 10)
        sut.open_run(0)
        sut.open_run(0)
    assert "already opened" in str(e.value)

def test_only_close_run():
    with pytest.raises(RuntimeError) as e:
        sut = Sequence("Foo", 10)
        sut.close_run(0)
    assert "already closed" in str(e.value)

def test_double_close_run():
    with pytest.raises(RuntimeError) as e:
        sut = Sequence("Foo", 1)
        sut.open_run(0)
        sut.close_run(1)
        sut.close_run(1)
    assert "already closed" in str(e.value)

def test_overlapping_runs():
    with pytest.raises(ValueError) as e:
        sut = Sequence("Foo", 19)
        sut.open_run(0)
        sut.close_run(15)

        sut.open_run(15)
        sut.close_run(19)
    assert "overlapping" in str(e.value)

def test_filter():
    # Arrange
    sequence = Sequence("Foo", 449)
    # Act
    sequence.open_run(0)
    sequence.close_run(150)
    sequence.open_run(300)
    sequence.close_run(449)
    sut = sequence.filter_runs(149)
    # Assert
    assert len(sut) == 1
    assert sut[0].start == 0
    assert sut[0].length == 150
    
def test_is_run_open():
    sut = Sequence("Foo", 150)
    sut.open_run(0)
    assert sut.is_run_open() == True
    sut.close_run(150)
    assert sut.is_run_open() == False

def test_negative_open():
    with pytest.raises(IndexError) as e:
        sut = Sequence("Foo", 10)
        sut.open_run(-10)
    assert "negative" in str(e.value)
    
def test_single_bucket():
    # Arrange
    sequence = Sequence("Foo", 1000)
    sequence.open_run(0)
    sequence.close_run(99)
    # Act
    # 1000/10 = 100 -> 10 buckets with lenght of 100
    sut = sequence.split_in_buckets(10)
    # Assert
    assert 0 in sut
    assert len(sut) == 1
    assert sut[0] == 100
    
def test_more_run_in_single_bucket():
    # Arrange
    sequence = Sequence("Foo", 1000)
    sequence.open_run(0)
    sequence.close_run(50)
    
    sequence.open_run(51)
    sequence.close_run(60)
    # Act
    # 1000/10 = 100 -> 10 buckets with lenght of 100
    sut = sequence.split_in_buckets(10)
    # Assert
    assert 0 in sut
    assert len(sut) == 1
    assert sut[0] == 61
    
def test_run_in_multiple_buckets():
    # Arrange
    sequence = Sequence("Foo", 1000)
    sequence.open_run(0)
    sequence.close_run(150)
    # Act
    # 1000/10 = 100 -> 10 buckets with lenght of 100
    sut = sequence.split_in_buckets(10)
    # Assert
    assert 0 in sut
    assert 1 in sut
    assert len(sut) == 2
    assert sut[0] == 100 # 0 -> 99
    assert sut[1] == 51 # 99 -> 151
    
def test_too_many_buckets():
    with pytest.raises(ValueError) as e:
        sut = Sequence("Foo", 1000)

        sut.open_run(0)
        sut.close_run(99)

        sut.open_run(100)
        sut.close_run(201)

        sut.split_in_buckets(1001)
    assert "at least" in str(e.value)
    
def test_run_end_to_end_divisible():
    sut = Sequence("Foo", 1000)

    sut.open_run(0)
    sut.close_run(1000)
    buckets = sut.split_in_buckets(10)
    for index in range(10):
        assert index in buckets
        assert buckets[index] == 100


def test_run_end_to_end_non_divisible():
    sut = Sequence("Foo", 1000)

    sut.open_run(0)
    sut.close_run(1000)
    buckets = sut.split_in_buckets(3)
    assert 0 in buckets
    assert 1 in buckets
    assert 2 in buckets    
    assert 3 in buckets
    
    assert len(buckets) == 4
    assert buckets[0] == 333
    assert buckets[1] == 333
    assert buckets[2] == 333
    assert buckets[3] == 2
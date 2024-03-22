from pathlib import Path
import pytest
from program.unknown_bases_stats import Sequence

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
        sut = Sequence("Foo", 10)
        sut.open_run(0)
        sut.close_run(1)
        sut.close_run(1)
    assert "already closed" in str(e.value)

def test_overlapping_runs():
    with pytest.raises(RuntimeError) as e:
        sut = Sequence("Foo", 10)
        sut.open_run(0)
        sut.close_run(15)

        sut.open_run(15)
        sut.close_run(19)
    assert "overlapping" in str(e.value)

def test_filter_runs_is_filtering():
    # Arrange
    sut = Sequence("Foo", 10)
    # Act
    sut.open_run(0)
    sut.close_run(150)
    sut.open_run(300)
    sut.close_run(449)
    runs = sut.filter_runs(149)
    # Assert
    assert len(runs) == 1
    assert runs[0].start == 0
    assert runs[0].length == 150
    

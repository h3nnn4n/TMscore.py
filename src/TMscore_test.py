import pytest
from TMscore import TMscore


def test_init():
    tmscore = TMscore("./TMscore")
    assert tmscore is not None


def test_init_no_args():
    tmscore = TMscore()
    assert tmscore is not None


def test_init_wrong_path():
    with pytest.raises(FileNotFoundError):
        TMscore('this path is bananas')


def test_call():
    tmscore = TMscore("./TMscore")
    tmscore("../proteins/1crn.pdb", "../proteins/best.pdb")
    data = tmscore.get_all()

    assert data['gdt_ha'] == (0.2772, (0.1304, 0.1957, 0.2609, 0.5217))
    assert data['gdt_ts'] == (0.4402, (0.1957, 0.2609, 0.5217, 0.7826))
    assert data['maxsub'] == 0.3286
    assert data['tm_score'] == 0.2875
    assert data['rmsd'] == 6.709


def test_get_rmsd():
    tmscore = TMscore("./TMscore")
    tmscore("../proteins/1crn.pdb", "../proteins/best.pdb")

    assert tmscore.get_rmsd() == 6.709


def test_get_tm_score():
    tmscore = TMscore("./TMscore")
    tmscore("../proteins/1crn.pdb", "../proteins/best.pdb")

    assert tmscore.get_tm_score() == 0.2875


def test_get_gdt_ts():
    tmscore = TMscore("./TMscore")
    tmscore("../proteins/1crn.pdb", "../proteins/best.pdb")

    assert tmscore.get_gdt_ts() == 0.4402


def test_get_gdt_ts_info():
    tmscore = TMscore("./TMscore")
    tmscore("../proteins/1crn.pdb", "../proteins/best.pdb")

    assert tmscore.get_gdt_ts_info() == (0.4402, (0.1957, 0.2609, 0.5217, 0.7826))


def test_get_gdt_ha():
    tmscore = TMscore("./TMscore")
    tmscore("../proteins/1crn.pdb", "../proteins/best.pdb")

    assert tmscore.get_gdt_ha() == 0.2772


def test_get_gdt_ha_info():
    tmscore = TMscore("./TMscore")
    tmscore("../proteins/1crn.pdb", "../proteins/best.pdb")

    assert tmscore.get_gdt_ha_info() == (0.2772, (0.1304, 0.1957, 0.2609, 0.5217))


def test_get_maxsub():
    tmscore = TMscore("./TMscore")
    tmscore("../proteins/1crn.pdb", "../proteins/best.pdb")

    assert tmscore.get_maxsub() == 0.3286

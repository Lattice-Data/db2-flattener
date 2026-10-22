"""GeneticModification alias labels for GEO column and guide filenames."""

import pytest

from db2_flattener.flatten.flattener import DB2Flattener


def make_flattener():
    return DB2Flattener.__new__(DB2Flattener)


@pytest.mark.parametrize(
    "aliases, expected",
    [
        (["lab:gm-screen-a"], "gm-screen-a"),
        ("lab:gm-screen-a", "gm-screen-a"),
        (["lab:foo/bar baz"], "foo_bar_baz"),
        (["no-prefix"], "no-prefix"),
        (["lab:keep.dots-and_underscores"], "keep.dots-and_underscores"),
    ],
)
def test_genetic_modification_alias_label(aliases, expected):
    assert (
        make_flattener()._genetic_modification_alias_label(
            aliases, gm_id="/genetic_modifications/x/"
        )
        == expected
    )


@pytest.mark.parametrize("aliases", [None, [], [""], ["  "]])
def test_genetic_modification_alias_label_rejects_empty(aliases):
    with pytest.raises(ValueError, match="exactly one alias"):
        make_flattener()._genetic_modification_alias_label(
            aliases, gm_id="/genetic_modifications/aaa/"
        )


def test_genetic_modification_alias_label_rejects_multiple():
    with pytest.raises(ValueError, match=r"/genetic_modifications/aaa/.*exactly one alias"):
        make_flattener()._genetic_modification_alias_label(
            ["lab:gm-a", "lab:gm-b"], gm_id="/genetic_modifications/aaa/"
        )


def test_validate_genetic_modification_aliases_ok():
    data = {
        "resolved_objects": {
            "GeneticModification": {
                "/genetic_modifications/aaa/": {
                    "@id": "/genetic_modifications/aaa/",
                    "aliases": ["lab:gm-a"],
                }
            }
        }
    }
    make_flattener()._validate_genetic_modification_aliases(data)


def test_validate_genetic_modification_aliases_raises():
    data = {
        "resolved_objects": {
            "GeneticModification": {
                "/genetic_modifications/aaa/": {
                    "@id": "/genetic_modifications/aaa/",
                    "aliases": ["lab:a", "lab:b"],
                }
            }
        }
    }
    with pytest.raises(ValueError, match="exactly one alias"):
        make_flattener()._validate_genetic_modification_aliases(data)

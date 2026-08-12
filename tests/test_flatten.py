from db2_flattener import flatten


def test_flatten_flat_dict():
    assert flatten({"a": 1, "b": 2}) == {"a": 1, "b": 2}


def test_flatten_nested_dict():
    assert flatten({"a": {"b": {"c": 1}}}) == {"a.b.c": 1}


def test_flatten_custom_separator():
    assert flatten({"a": {"b": 1}}, sep="/") == {"a/b": 1}


def test_flatten_empty():
    assert flatten({}) == {}

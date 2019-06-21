# Example reduction scripts

The scripts are written as simple unit tests to move the discussion in that direction.
I chose to use python's `unittest` for simplicity.

The tests should come in pairs, one for the old API and a corresponding one for the mock API.
Both should give the same result.

The goal is to be able to import from the real new API and run the same scripts.
More details on that will come as the new API is made available.

To run the tests, simply do:

```
python test_mock_api.py
```
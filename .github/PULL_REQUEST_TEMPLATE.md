## What this changes

<!-- One or two sentences. Link the issue it closes, if there is one. -->

## Why

<!-- What was wrong, or what this makes possible. -->

## How it was tested

<!--
Be specific. "Ran the deconstruction over 200 CSD structures, 3 previously
failing ones now resolve" tells a reviewer far more than "should work".
Name the structures you tried it on.
-->

## Checklist

- [ ] `pytest tests/` passes
- [ ] A test covers the change, with the structure that reproduced it added to
      `tests/test_data/` if this fixes a bug
- [ ] `CHANGELOG.md` updated under the next release heading
- [ ] Docstrings added or updated for any new public function
- [ ] Docs still build without warnings (`sphinx-build -b html docs/source docs/build/html`)

# Contributing to mofstructure

Thanks for taking the time to contribute. Bug reports, structure files that
break the deconstruction, and pull requests are all welcome.

## Reporting a bug

Most problems in `mofstructure` are specific to one framework rather than
general, so **please attach the structure file**. A CIF that reproduces the
problem is worth more than any description of it, and without one a report
usually cannot be acted on.

Open an issue at
[github.com/bafgreat/mofstructure/issues](https://github.com/bafgreat/mofstructure/issues)
and include:

- the structure file, or a link to it in the CSD or elsewhere
- the command or Python call you ran
- what you expected and what happened
- output of `pip show mofstructure` and your Python version
- `java -version` if the problem involves topology

## Suggesting a feature

Open an issue describing the analysis you need and, where possible, the
chemistry behind it. Rules for deconstructing an unusual coordination
environment are much easier to implement with a concrete example to test
against.

## Development setup

```bash
git clone https://github.com/bafgreat/mofstructure.git
cd mofstructure
pip install -e .
```

Run the tests before opening a pull request:

```bash
pytest tests/
```

Topology requires a JRE. Everything else runs without one.

## Pull requests

- Branch from `main` and keep one logical change per pull request.
- Add a test when you fix a bug. `tests/test_data/` holds small structures for
  this; add the CIF that reproduced the problem.
- Match the surrounding style. Docstrings use `'''` with `**parameters:**` and
  `**returns:**` sections.
- Update `CHANGELOG.md` under a heading for the next release.
- Say in the description what you tested it on. "Ran the full deconstruction
  over 200 CSD structures" tells a reviewer far more than "should work".

### Things worth knowing before you change the analysis code

- **SMILES are toolkit specific.** OpenBabel, RDKit and CACTVS each produce a
  different canonical SMILES for the same molecule. Never compare SMILES across
  toolkits, and never use them as dictionary keys without normalising both
  sides through `mofdeconstructor.name_lookup_keys()`.
- **Ligands come out of deconstruction with open valences.** A fragment cut
  from a framework is not the neutral parent molecule and will not match a
  reference database until its valences are saturated.
- **zeo++ can abort the process.** It calls `abort()` on structures whose
  Voronoi decomposition fails an internal check, which no `try`/`except` can
  catch. This is why `porosity.zeo_calculation` runs in a child interpreter.
  Do not move it back in process.
- **The CLI skips structures it has already done.** Delete the output folder
  when testing a change, or you will be looking at cached results.

## Code of conduct

Be civil and constructive. Discussion should stay about the code and the
chemistry.

## Licence

Contributions are accepted under the MIT Licence, the same terms as the
project. See [LICENSE](LICENSE).

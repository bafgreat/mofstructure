'''
mofstructure: deconstruction, porosity and topology of porous frameworks.

The package takes a periodic structure and answers three kinds of question
about it. `structure.MOFstructure` is the single entry point for most work:
it removes unbound guests, deconstructs a framework into its building units,
and reports cheminformatic identifiers for each. `porosity` measures the pore
geometry through zeo++. `topology` names the underlying net, working out for
itself whether it has been handed a MOF, a COF or a zeolite.

The same analyses are available from the command line; see the console
scripts in `mofstructure.scripts`.
'''

# Makefile before PR
#

.PHONY: checks

checks: flake spellcheck

flake:
	flake8

spellcheck:
	codespell buildtest/ docs/ examples/

test:
	pytest -v
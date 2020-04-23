.PHONY: use_cases.pdf all clean

all: use_cases.pdf

use_cases.pdf: use_cases.tex library.bib
	latexmk -pdf use_cases.tex

clean:
	latexmk -pdf -CA

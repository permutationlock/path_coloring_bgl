LATEX = pdflatex
NAME = path_list_coloring_examples
    # Base for all filenames

all: $(NAME).pdf
.PHONY: all

$(NAME).pdf: $(NAME).tex
	$(LATEX) $(NAME).tex
	$(LATEX) $(NAME).tex

clean:
	rm -f $(NAME).pdf $(NAME).aux $(NAME).log $(NAME).out $(NAME).pre
.PHONY: clean


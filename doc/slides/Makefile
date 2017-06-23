LATEX = pdflatex
NAME = presentation
    # Base for all filenames

all: $(NAME).pdf
.PHONY: all

$(NAME).pdf: $(NAME).tex
	$(LATEX) $(NAME).tex
	$(LATEX) $(NAME).tex

clean:
	rm -f $(NAME).pdf $(NAME).aux $(NAME).log $(NAME).out $(NAME).pre \
	$(NAME).nav $(NAME).snm $(NAME).toc
.PHONY: clean


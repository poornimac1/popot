PLOT_SCRIPTS = $(wildcard *.plot)
GIF_FROM_PLOT=$(patsubst %.plot, %.gif, $(PLOT_SCRIPTS))
PDF_FROM_PLOT = $(patsubst %.plot, %.pdf, $(PLOT_SCRIPTS))
MOVIE=movie.avi
ANIMATION=movie.gif

%.gif : %.plot
	@gnuplot $<
	@echo $@ "generated."

%.fig : %.plot
	@gnuplot $<
	@echo $@ "generated."

%.pdf : %.fig
	@fig2dev -L pdf $< $@
	@echo $@ "generated."

clean:
	rm -f *.data *.plot *.gif *~ *.pdf 

$(MOVIE):$(GIF_FROM_PLOT)
	ffmpeg -i %05d.gif -b 1M $(MOVIE)

$(ANIMATION):$(GIF_FROM_PLOT)
	convert -resize 50% -delay 5 -loop 0 0*.gif $(ANIMATION)

all:$(ANIMATION)

PLOT_SCRIPTS = $(wildcard *.plot)
GIF_FROM_PLOT=$(patsubst %.plot, %.gif, $(PLOT_SCRIPTS))
CONNEC_SCRIPTS = $(wildcard *.dot)
GIF_FROM_CONNEC=$(patsubst %.dot, %.gif_connec, $(CONNEC_SCRIPTS))
MOVIE=movie.avi
ANIMATION=movie.gif
ANIMATION_CONNECTIONS=movie_connections.gif

%.gif : %.plot
	@gnuplot $<
	@echo $@ "generated."

%.gif_connec : %.dot
	@neato -Tgif $< > $@
	@echo $@ "generated."

clean:
	rm -f *.data *.plot *.gif *~ *.pdf 

$(MOVIE):$(GIF_FROM_PLOT)
	ffmpeg -i %05d.gif -b 1M $(MOVIE)

$(ANIMATION):$(GIF_FROM_PLOT)
	convert -resize 50% -delay 5 -loop 0 0*.gif $(ANIMATION)

$(ANIMATION_CONNECTIONS):$(GIF_FROM_CONNEC)
	convert -resize 25% -delay 5 -loop 0 0*.gif_connec $(ANIMATION_CONNECTIONS)

all:$(ANIMATION) $(ANIMATION_CONNECTIONS)

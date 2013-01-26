G77 = g77

all: gerg2004

gerg2004: gerg2004.f
	$(G77) $< -o $@

clean:
	rm -rf gerg2004
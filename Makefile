all:
	make -C src all

install:
	R CMD INSTALL .

clean:
	make -C src clean

build:
	R CMD build .

check: build
	R CMD check --no-manual `ls -1tr geiger*gz | tail -n1`
	@rm -f `ls -1tr geiger*gz | tail -n1`
	@rm -rf geiger.Rcheck

.PHONY: all install clean check

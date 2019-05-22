# GOS_NMF
## Target
current: target
-include target.mk

## Makestuff setup
msrepo = https://github.com/dushoff
ms = makestuff
-include $(ms)/os.mk

Ignore += $(ms)
Makefile: $(ms)/Makefile
$(ms)/Makefile:
	git clone $(msrepo)/$(ms)
	ls $@

######################################################################

Sources += Makefile README.md

# README.html: README.md
Ignore += README.html

Sources += $(wildcard *.R)

######################################################################

### Makestuff rules

-include $(ms)/pandoc.mk
-include $(ms)/git.mk
-include $(ms)/visual.mk


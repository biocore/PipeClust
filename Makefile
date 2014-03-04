#-------------------------
# PipeClust Makefile
#
# Author:	Jose Antonio Navas Molina
# Contact:	josenavasmolina@gmail.com
# Date: 	2014-03-03
#
#-------------------------

TARGET	= PipeClust

SRCDIR	= src
INCLDIR	= include
OBJDIR	= obj
BINDIR	= bin

SOURCES		:= $(wildcard $(SRCDIR)/*.c)
INCLUDES	:= $(wildcard $(INCLDIR)/*.h)
OBJECTS		:= $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

CC		= mpicc
CFLAGS	= -Wall -O3 -msse2 -c -I ./${INCLDIR}/
# CFLAGS	= -Wall -g -pg -c -I ./${INCLDIR}/

LINKER	= mpicc
LFLAGS	= -Wall -O3 -msse2 -lm
# LFLAGS	= -Wall -g -pg -lm

MKDIR	= mkdir -p

$(BINDIR)/$(TARGET): $(OBJECTS) ${BINDIR}
	$(CC) -o $@ $(LFLAGS) $(OBJECTS)

${BINDIR}:
	${MKDIR} ${BINDIR}

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c ${OBJDIR}
	$(CC) $(CFLAGS) $< -o $@

${OBJDIR}:
	${MKDIR} ${OBJDIR}

.PHONY: clean

clean:
	rm -rf $(OBJDIR) $(BINDIR)

.PHONY: clean_obj

clean_obj:
	rm -rf $(OBJDIR)
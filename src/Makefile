CPPSRCS= main.cc cmdlineopts.cc driver.cc phaser.cc
CSRCS= # admutils.c mcio.c # files for Nick's I/O
OBJS= $(patsubst %.cc,%.o,$(CPPSRCS)) $(patsubst %.c,%.o,$(CSRCS))
EXEC= hapi-ur

GPP = g++
GCC = gcc
#DEFINES= -DPROFILE
#DEFINES= -DIBD
#DEFINES= -DSTATS
DEFINES= 
INCLUDES= -I../include -I../lib -I../sparsehash/include
#CFLAGS = -g -Wall $(DEFINES) $(INCLUDES)
# optimized; remove asserts
#CFLAGS = -O2 -Wall -DNDEBUG $(DEFINES) $(INCLUDES)
#CFLAGS = -O2 -Wall $(DEFINES) $(INCLUDES)
CFLAGS = -O2 -Wall -g $(DEFINES) $(INCLUDES)
# profiling:
#CFLAGS = -pg -O2 -Wall $(DEFINES) $(INCLUDES)

LIB_PATH= -L../lib
# librt for clock_gettime()
LIBS = -lrt -lgenetio $(LIB_PATH)

# dependency variables / commands
DEPDIR = .deps
df = $(DEPDIR)/$(*F)

all: $(EXEC)

$(EXEC): $(OBJS) $(HEADERS)
	$(GPP) -o $(EXEC) $(OBJS) $(CFLAGS) $(LIBS)

vote-phase: vote-phase.cc
	$(GPP) -o $@ -Wall -O2 $^

# This way of building dependencies (per-file) described at
# http://make.paulandlesley.org/autodep.html

.c.o:
	@mkdir -p $(DEPDIR)
	$(GCC) -MMD $(CFLAGS) -o $@ -c $<
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d

.cc.o:
	@mkdir -p $(DEPDIR)
	$(GPP) -MMD $(CFLAGS) -o $@ -c $<
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d


# include the .P dependency files, but don't warn if they don't exist (the -)
-include $(CPPSRCS:%.cc=$(DEPDIR)/%.P)
-include $(CSRCS:%.c=$(DEPDIR)/%.P)
# The following applies if we don't use a dependency directory:
#-include $(SRCS:.cc=.P)

tags: $(SRCS) *.h
	ctags --language-force=c++ --extra=+q --fields=+i --excmd=n *.cc *.h

clean:
	rm -f $(EXEC) $(OBJS) vote-phase

clean-deps:
	rm -f $(DEPDIR)/*.P

CPPSRCS= marker.cc util.cc hapi-ur-util.cc personio.cc personbits.cc
CSRCS= # admutils.c mcio.c # files for Nick's I/O
OBJS= $(patsubst %.cc,%.o,$(CPPSRCS)) $(patsubst %.c,%.o,$(CSRCS))
LIB=  ../libgenetio.a

GPP = g++
GCC = gcc
AR  = ar
DEFINES= 
INCLUDES= -I../../include
#CFLAGS = -g -Wall $(DEFINES)
# optimized; remove asserts
#CFLAGS = -O2 -Wall -DNDEBUG $(DEFINES)
CFLAGS = -O2 -Wall $(DEFINES) $(INCLUDES)
# profiling:
#CFLAGS = -pg -O2 -Wall $(DEFINES)

# dependency variables / commands
DEPDIR = .deps
df = $(DEPDIR)/$(*F)

all: $(LIB)

$(LIB): $(OBJS) $(HEADERS)
	$(AR) -cvq $(LIB) $(OBJS)

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
	ctags --language-force=c++ --extra=+q --fields=+i --excmd=n *.c *.cc *.h

clean:
	rm -f $(LIB) $(OBJS)

clean-deps:
	rm -f $(DEPDIR)/*.P

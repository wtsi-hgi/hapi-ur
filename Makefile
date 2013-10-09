SUBDIRS = lib/genetio src
ECHO = echo

all : lib/libgenetio.a src/hapi-ur src/vote-phase

lib/libgenetio.a : force_look
	$(ECHO) looking into lib/genetio : $(MAKE) $(MFLAGS)
	cd lib/genetio; $(MAKE) $(MFLAGS)

src/hapi-ur src/vote-phtase : force_look
	$(ECHO) looking into src : $(MAKE) $(MFLAGS)
	cd src; $(MAKE) $(MFLAGS)

clean : 
	-for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean ); done

force_look : 
	true

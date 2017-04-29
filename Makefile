CXX=g++
RM=rm -f
CPPFLAGS=-g -I. -std=c++11 -O3
LDFLAGS=-g
LDLIBS=-static-libgcc -static-libstdc++

SRCS:=$(shell find . -name "*.cpp")
OBJS=$(subst .cpp,.o,$(SRCS))

all: srp2

srp2: $(OBJS)
	$(CXX)  $(LDFLAGS) -o srp2 $(OBJS) $(LDLIBS) 

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~ .depend

include .depend
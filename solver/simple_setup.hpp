#ifndef _SIMPLE_SETUP_H_
#define _SIMPLE_SETUP_H_

template <typename Lattice> class NoSetup {
public:
  NoSetup(Lattice&) {}

  void setup() const {}
};


#endif /* _SIMPLE_SETUP_H_ */


#ifndef MEDYAN_CAMKII_Species_h
#define MEDYAN_CAMKII_Species_h

class SpeciesCaMKIIer : public SpeciesBound {

public:
	/// Default constructor
	SpeciesCaMKIIer() : SpeciesBound() {}

	/// The main constructor
	/// @param name - Example, "G-Actin" or "Arp2/3"
	/// @param n - copy number
	SpeciesCaMKIIer (const string &name, species_copy_t n=0, species_copy_t ulim=1)
			:  SpeciesBound(name, n, ulim) {};

	/// Copy constructor
	SpeciesCaMKIIer(const SpeciesCaMKIIer &rhs)  : SpeciesBound(rhs) {}

	/// Move constructor
	SpeciesCaMKIIer (SpeciesCaMKIIer &&rhs) noexcept : SpeciesBound(move(rhs)){
	}

	/// Regular Assignment
	SpeciesCaMKIIer& operator=(const SpeciesCaMKIIer& rhs)  {
		SpeciesBound::operator=(rhs);
		return *this;
	}

	/// Move assignment
	SpeciesCaMKIIer& operator=(SpeciesCaMKIIer&& rhs)
	{
		SpeciesBound::operator=(move(rhs));
		return *this;
	}

	virtual SpeciesCaMKIIer* clone() {
		return new SpeciesCaMKIIer(*this);
	}

	/// Return the full name of this Species in a string format
	virtual string getFullName() const {return getName() + "{CaMKIIer}";}

	/// Default destructor
	~SpeciesCaMKIIer () noexcept {};
};


class SpeciesCaMKIICylinder : public SpeciesBound {

public:
	/// Default constructor
	SpeciesCaMKIICylinder() : SpeciesBound() {}

	/// The main constructor
	/// @param name - Example, "G-Actin" or "Arp2/3"
	/// @param n - copy number
	SpeciesCaMKIICylinder (const string &name, species_copy_t n=0, species_copy_t ulim=1)
			:  SpeciesBound(name, n, ulim) {};

	/// Copy constructor
	SpeciesCaMKIICylinder(const SpeciesCaMKIICylinder &rhs)  : SpeciesBound(rhs) {}

	/// Move constructor
	SpeciesCaMKIICylinder (SpeciesCaMKIICylinder &&rhs) noexcept : SpeciesBound(move(rhs)){
	}

	/// Regular Assignment
	SpeciesCaMKIICylinder& operator=(const SpeciesCaMKIICylinder& rhs)  {
		SpeciesBound::operator=(rhs);
		return *this;
	}

	/// Move assignment
	SpeciesCaMKIICylinder& operator=(SpeciesCaMKIICylinder&& rhs)
	{
		SpeciesBound::operator=(move(rhs));
		return *this;
	}

	virtual SpeciesCaMKIICylinder* clone() {
		return new SpeciesCaMKIICylinder(*this);
	}

	/// Return the full name of this Species in a string format
	virtual string getFullName() const {return getName() + "{CaMKIICylinderSpecies}";}

	/// Default destructor
	~SpeciesCaMKIICylinder () noexcept {};
};



#endif //MEDYAN_CAMKII_Species_h

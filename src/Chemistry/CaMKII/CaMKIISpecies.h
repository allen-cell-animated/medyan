
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

	/// Return the full name of this Species in a string format (e.g. "Arp2/3{CaMKIIer}"
	virtual string getFullName() const {return getName() + "{CaMKIIer}";}

	/// Default destructor
	~SpeciesCaMKIIer () noexcept {};
};


class SpeciesCaMKIIDummyCylinder : public SpeciesBound {

public:
	/// Default constructor
	SpeciesCaMKIIDummyCylinder() : SpeciesBound() {}

	/// The main constructor
	/// @param name - Example, "G-Actin" or "Arp2/3"
	/// @param n - copy number
	SpeciesCaMKIIDummyCylinder (const string &name, species_copy_t n=0, species_copy_t ulim=1)
			:  SpeciesBound(name, n, ulim) {};

	/// Copy constructor
	SpeciesCaMKIIDummyCylinder(const SpeciesCaMKIIDummyCylinder &rhs)  : SpeciesBound(rhs) {}

	/// Move constructor
	SpeciesCaMKIIDummyCylinder (SpeciesCaMKIIDummyCylinder &&rhs) noexcept : SpeciesBound(move(rhs)){
	}

	/// Regular Assignment
	SpeciesCaMKIIDummyCylinder& operator=(const SpeciesCaMKIIDummyCylinder& rhs)  {
		SpeciesBound::operator=(rhs);
		return *this;
	}

	/// Move assignment
	SpeciesCaMKIIDummyCylinder& operator=(SpeciesCaMKIIDummyCylinder&& rhs)
	{
		SpeciesBound::operator=(move(rhs));
		return *this;
	}

	virtual SpeciesCaMKIIDummyCylinder* clone() {
		return new SpeciesCaMKIIDummyCylinder(*this);
	}

	/// Return the full name of this Species in a string format (e.g. "Arp2/3{CaMKIIer}"
	virtual string getFullName() const {return getName() + "{CaMKIIDummyCylinder}";}

	/// Default destructor
	~SpeciesCaMKIIDummyCylinder () noexcept {};
};



#endif //MEDYAN_CAMKII_Species_h

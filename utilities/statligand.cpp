#include <iostream>
#include <iomanip>
#include <array>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <limits>
#include <cassert>
using namespace std;

class atom
{
private:
	static const size_t n = 31;
	static const array<string, n> ad_strings;
	static const array<float, n> ad_covalent_radii;
	static const array<float, n> ad_atomic_weights;
	static const array<size_t, n> ad_to_xs;
public:
	size_t serial;
	array<float, 3> coord;
	size_t ad;
	size_t xs;

	// Constructs an atom from an ATOM/HETATM line in PDBQT format.
	explicit atom(const string& line);

	/// Returns true if the AutoDock4 atom type is not supported.
	bool ad_unsupported() const;

	/// Returns true if the atom is a nonpolar hydrogen atom.
	bool is_nonpolar_hydrogen() const;

	/// Returns true if the atom is a polar hydrogen atom.
	bool is_polar_hydrogen() const;

	/// Returns the covalent radius of current AutoDock4 atom type.
	float covalent_radius() const;

	/// Returns the atomic weight of current AutoDock4 atom type.
	float atomic_weight() const;

	/// Returns true if the current atom is covalently bonded to a given atom.
	bool has_covalent_bond(const atom& a) const;
};

/// AutoDock4 atom type names.
const array<string, atom::n> atom::ad_strings =
{
	"H" , //  0
	"HD", //  1
	"C" , //  2
	"A" , //  3
	"N" , //  4
	"NA", //  5
	"OA", //  6
	"S" , //  7
	"SA", //  8
	"Se", //  9
	"P" , // 10
	"F" , // 11
	"Cl", // 12
	"Br", // 13
	"I" , // 14
	"Zn", // 15
	"Fe", // 16
	"Mg", // 17
	"Ca", // 18
	"Mn", // 19
	"Cu", // 20
	"Na", // 21
	"K" , // 22
	"Hg", // 23
	"Ni", // 24
	"Co", // 25
	"Cd", // 26
	"As", // 27
	"Sr", // 28
	"U" , // 29
	"Cs", // 30
};

/// AutoDock4 covalent radii, factorized by 1.1 for extra allowance.
const array<float, atom::n> atom::ad_covalent_radii =
{
	0.407f, //  0 = H , 0.407 = 1.1 * 0.37
	0.407f, //  1 = HD, 0.407 = 1.1 * 0.37
	0.847f, //  2 = C , 0.847 = 1.1 * 0.77
	0.847f, //  3 = A , 0.847 = 1.1 * 0.77
	0.825f, //  4 = N , 0.825 = 1.1 * 0.75
	0.825f, //  5 = NA, 0.825 = 1.1 * 0.75
	0.803f, //  6 = OA, 0.803 = 1.1 * 0.73
	1.122f, //  7 = S , 1.122 = 1.1 * 1.02
	1.122f, //  8 = SA, 1.122 = 1.1 * 1.02
	1.276f, //  9 = Se, 1.276 = 1.1 * 1.16
	1.166f, // 10 = P , 1.166 = 1.1 * 1.06
	0.781f, // 11 = F , 0.781 = 1.1 * 0.71
	1.089f, // 12 = Cl, 1.089 = 1.1 * 0.99
	1.254f, // 13 = Br, 1.254 = 1.1 * 1.14
	1.463f, // 14 = I , 1.463 = 1.1 * 1.33
	1.441f, // 15 = Zn, 1.441 = 1.1 * 1.31
	1.375f, // 16 = Fe, 1.375 = 1.1 * 1.25
	1.430f, // 17 = Mg, 1.430 = 1.1 * 1.30
	1.914f, // 18 = Ca, 1.914 = 1.1 * 1.74
	1.529f, // 19 = Mn, 1.529 = 1.1 * 1.39
	1.518f, // 20 = Cu, 1.518 = 1.1 * 1.38
	1.694f, // 21 = Na, 1.694 = 1.1 * 1.54
	2.156f, // 22 = K , 2.156 = 1.1 * 1.96
	1.639f, // 23 = Hg, 1.639 = 1.1 * 1.49
	1.331f, // 24 = Ni, 1.331 = 1.1 * 1.21
	1.386f, // 25 = Co, 1.386 = 1.1 * 1.26
	1.628f, // 26 = Cd, 1.628 = 1.1 * 1.48
	1.309f, // 27 = As, 1.309 = 1.1 * 1.19
	2.112f, // 28 = Sr, 2.112 = 1.1 * 1.92
	2.156f, // 29 = U , 2.156 = 1.1 * 1.96
	2.475f, // 30 = Cs, 2.475 = 1.1 * 2.25
};

/// AutoDock4 atomic weights.
const array<float, atom::n> atom::ad_atomic_weights =
{
	  1.008,//  0 = H
	  1.008,//  1 = HD
	 12.01, //  2 = C
	 12.01, //  3 = A
	 14.01, //  4 = N
	 14.01, //  5 = NA
	 16.00, //  6 = OA
	 32.07, //  7 = S
	 32.07, //  8 = SA
	 78.96, //  9 = Se
	 30.97, // 10 = P
	 19.00, // 11 = F
	 35.45, // 12 = Cl
	 79.90, // 13 = Br
	126.90, // 14 = I
	 65.38, // 15 = Zn
	 55.85, // 16 = Fe
	 24.31, // 17 = Mg
	 40.08, // 18 = Ca
	 54.94, // 19 = Mn
	 63.55, // 20 = Cu
	 22.99, // 21 = Na
	 39.10, // 22 = K
	200.59, // 23 = Hg
	 58.69, // 24 = Ni
	 58.93, // 25 = Co
	112.41, // 26 = Cd
	 74.92, // 27 = As
	 87.62, // 28 = Sr
	238.03, // 29 = U
	132.91, // 30 = Cs
};

/// Mapping from AutoDock4 atom type to XScore atom type.
const array<size_t, atom::n> atom::ad_to_xs =
{
	 n, //  0 = H  -> dummy
	 n, //  1 = HD -> dummy
	 0, //  2 = C  -> C_H   =  0, Carbon, hydrophobic, not bonded to a hetero atom.
	 0, //  3 = A  -> C_H   =  0, Carbon, hydrophobic, not bonded to a hetero atom.
	 2, //  4 = N  -> N_P   =  2, Nitrogen, neither hydrogen bond donor nor acceptor.
	 4, //  5 = NA -> N_A   =  4, Nitrogen, hydrogen bond acceptor.
	 6, //  6 = OA -> O_A   =  6, Oxygen, hydrogen bond acceptor.
	 8, //  7 = S  -> S_P   =  8, Sulfur or Selenium.
	 8, //  8 = SA -> S_P   =  8, Sulfur or Selenium.
	 8, //  9 = Se -> S_P   =  8, Sulfur or Selenium.
	 9, // 10 = P  -> P_P   =  9, Phosphorus.
	10, // 11 = F  -> F_H   = 10, Fluorine, hydrophobic.
	11, // 12 = Cl -> Cl_H  = 11, Chlorine, hydrophobic.
	12, // 13 = Br -> Br_H  = 12, Bromine, hydrophobic.
	13, // 14 = I  -> I_H   = 13, Iodine, hydrophobic.
	14, // 15 = Zn -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 16 = Fe -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 17 = Mg -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 18 = Ca -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 19 = Mn -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 20 = Cu -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 21 = Na -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 22 = K  -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 23 = Hg -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 24 = Ni -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 25 = Co -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 26 = Cd -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 27 = As -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 28 = Sr -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 29 = U  -> Met_D = 14, Metal, hydrogen bond donor.
	14, // 30 = Cs -> Met_D = 14, Metal, hydrogen bond donor.
};

atom::atom(const string& line) :
	ad(find(ad_strings.cbegin(), ad_strings.cend(), line.substr(77, isspace(line[78]) ? 1 : 2)) - ad_strings.cbegin()),
	xs(ad_to_xs[ad])
{
	coord[0] = stof(line.substr(30, 8));
	coord[1] = stof(line.substr(38, 8));
	coord[2] = stof(line.substr(46, 8));
}

/// Returns true if the AutoDock4 atom type is not supported.
bool atom::ad_unsupported() const
{
	return ad == n;
}

/// Returns true if the atom is a nonpolar hydrogen atom.
bool atom::is_nonpolar_hydrogen() const
{
	return ad == 0;
}

/// Returns true if the atom is a polar hydrogen atom.
bool atom::is_polar_hydrogen() const
{
	return ad == 1;
}

/// Returns the covalent radius of current AutoDock4 atom type.
float atom::covalent_radius() const
{
	return ad_covalent_radii[ad];
}

/// Returns the atomic weight of current AutoDock4 atom type.
float atom::atomic_weight() const
{
	return ad_atomic_weights[ad];
}

bool atom::has_covalent_bond(const atom& a) const
{
	const float d0 = coord[0] - a.coord[0];
	const float d1 = coord[1] - a.coord[1];
	const float d2 = coord[2] - a.coord[2];
	const float s = covalent_radius() + a.covalent_radius();
	return d0 * d0 + d1 * d1 + d2 * d2 < s * s;
}

//! Returns true if the XScore atom type is a hydrogen bond donor.
inline bool is_hbdonor(const size_t t)
{
	return t == 14;
}

//! Returns true if the XScore atom type is a hydrogen bond acceptor.
inline bool is_hbacceptor(const size_t t)
{
	return t ==  4 || t ==  5 || t ==  6 || t ==  7;
}

class ligand : public vector<atom>
{
public:
	/// Load current ligand from an ifstream
	explicit ligand(istream& ifs);

	/// Ligand properties.
	size_t num_hydrogens;
	size_t num_hydrogen_bond_donors;
	size_t num_hydrogen_bond_acceptors;
	size_t num_active_torsions;
	size_t num_inactive_torsions;
	float molecular_weight;
	vector<float> mn;
	vector<float> mx;
};

/// Represents a ROOT or a BRANCH in PDBQT structure.
class frame
{
public:
	size_t parent; ///< Frame array index pointing to the parent of current frame. For ROOT frame, this field is not used.
	size_t rotorYidx; ///< Index pointing to the current frame atom which forms a rotatable bond with the rotorX atom of parent frame.

	/// Constructs an active frame, and relates it to its parent frame.
	explicit frame(const size_t parent, const size_t rotorYidx) : parent(parent), rotorYidx(rotorYidx) {}
};

ligand::ligand(istream& ifs) : num_hydrogens(0), num_hydrogen_bond_donors(0), num_hydrogen_bond_acceptors(0), num_active_torsions(0), num_inactive_torsions(0), molecular_weight(0), mn(3, numeric_limits<float>::max()), mx(3, numeric_limits<float>::lowest())
{
	// Initialize necessary variables for constructing a ligand.
	vector<frame> frames; ///< ROOT and BRANCH frames.
	frames.reserve(30); // A ligand typically consists of <= 30 frames.
	frames.emplace_back(0, 0); // ROOT is also treated as a frame. The parent of ROOT frame is dummy.
	reserve(100); // A ligand typically consists of <= 100 heavy atoms.

	// Initialize helper variables for parsing.
	size_t current = 0; // Index of current frame, initialized to ROOT frame.
	frame* f = &frames.front(); // Pointer to the current frame.

	// Parse the ligand line by line.
	for (string line; getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			// Parse the line.
			atom a(line);

			// Skip unsupported atom types.
			if (a.ad_unsupported()) continue;

			// Count statistics.
			if (a.is_nonpolar_hydrogen())
			{
				++num_hydrogens;
			}
			else if (a.is_polar_hydrogen())
			{
				++num_hydrogens;
				++num_hydrogen_bond_donors;
			}
			else // Current atom is a heavy atom.
			{
				num_hydrogen_bond_donors += is_hbdonor(a.xs);
				num_hydrogen_bond_acceptors += is_hbacceptor(a.xs);
				// Save the heavy atom.
				push_back(move(a));
			}
			molecular_weight += a.atomic_weight();
			for (size_t i = 0; i < 3; ++i)
			{
				mn[i] = min<float>(mn[i], a.coord[i]);
				mx[i] = max<float>(mx[i], a.coord[i]);
			}
		}
		else if (record == "BRANCH")
		{
			// Insert a new frame whose parent is the current frame.
			frames.push_back(frame(current, size()));

			// Now the current frame is the newly inserted BRANCH frame.
			current = frames.size() - 1;

			// Update the pointer to the current frame.
			f = &frames[current];
		}
		else if (record == "ENDBRA")
		{
			// If the current frame consists of rotor Y and a few hydrogens only, e.g. -OH, -NH2 or -CH3,
			// the torsion of this frame will have no effect on scoring and is thus redundant.
			if (current + 1 == frames.size() && f->rotorYidx + 1 == size())
			{
				++num_inactive_torsions;
			}
			else
			{
				++num_active_torsions;
			}

			// Now the parent of the following frame is the parent of current frame.
			current = f->parent;

			// Update the pointer to the current frame.
			f = &frames[current];
		}
		else if (record == "TORSDO") break;
	}
	assert(1 + num_active_torsions + num_inactive_torsions == frames.size());
	assert(num_hydrogen_bond_acceptors <= size());
	assert(num_hydrogen_bond_donors + num_hydrogen_bond_acceptors <= num_hydrogens + size());
}

int main(int argc, char* argv[])
{
	ligand lig(cin);
	if (lig.empty())
	{
		cout << "H,HA,HBD,HBA,NAT,NIT,MWT,corner0_x,corner0_y,corner0_z,corner1_x,corner1_y,corner1_z" << endl;
	}
	else
	{
		cout.setf(ios::fixed, ios::floatfield);
		cout << lig.num_hydrogens + lig.size() << ',' << lig.size() << ',' << lig.num_hydrogen_bond_donors << ',' << lig.num_hydrogen_bond_acceptors << ',' << lig.num_active_torsions << ',' << lig.num_inactive_torsions << ',' << setprecision(3) << lig.molecular_weight << ',' << lig.mn[0] << ',' << lig.mn[1] << ',' << lig.mn[2] << ',' << lig.mx[0] << ',' << lig.mx[1] << ',' << lig.mx[2] << endl;
	}
}

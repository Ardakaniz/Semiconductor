#include <iostream>
#include <fstream>
#include <vector>

#include "Quantity.hpp"

void iter(Quantity& n, Quantity& p, const Quantity& Cd, const Quantity& Ca, std::vector<double>& electric_field, unsigned int t) {
	const double K = 1.0; //* == \mu * e / \epsilon

	double qty_integral = n(t, 0) - p(t, 0) + (Ca(t, 0) - Cd(t, 0));
	for (int j = 0; j < static_cast<int>(n.get_pos_count()); ++j) {
		const double dn = 0.5 * (n(t, j + 1) - n(t, j - 1));
		const double dp = 0.5 * (p(t, j + 1) - p(t, j - 1));
		const double ddn = 0;//(n(t, j + 1) + n(t, j - 1) - 2 * n(t, j)) / (n.get_pos_step() * n.get_pos_step());
		const double ddp = 0;//(p(t, j + 1) + p(t, j - 1) - 2 * p(t, j)) / (p.get_pos_step() * p.get_pos_step());

		const double E = 0.5 * (p(t, 0) + p(t, j) - n(t, 0) - n(t, j) + Cd(t, 0) + Cd(t, j) - Ca(t, 0) - Ca(t, j)) + qty_integral;
		const double total_charges = p(t, j) - n(t, j) + (Cd(t, j) - Ca(t, j));
		
		electric_field.push_back(E * n.get_pos_step());
		n(t + 1, j) = n(t, j) + n.get_time_step() * (K * (dn * E + total_charges * n(t, j)) + ddn);
		p(t + 1, j) = p(t, j) - p.get_time_step() * (K * (dp * E + total_charges * p(t, j)) - ddp);

		qty_integral += total_charges;
	}
}

int main() {
	const double time_end = 1;
	const unsigned int time_count = 100000;
	const double time_step = time_end / static_cast<double>(time_count);
	const double pos_end = 1.0;
	const unsigned int pos_count = 100; // 10000
	const double pos_step = pos_end / static_cast<double>(pos_count);


	Quantity n{ time_step, time_count, pos_step, pos_count, true };
	Quantity p{ time_step, time_count, pos_step, pos_count, true };
	Quantity Cd{ time_step, time_count, pos_step, pos_count };
	Quantity Ca{ time_step, time_count, pos_step, pos_count };

	const double n_i = 0.001; // initial charge density
	
	for (unsigned int j = 0; j < pos_count; ++j) {
		if (j < pos_count / 2 - 20) {
			p(0, j) = n_i;
		}
		else if (j >= pos_count / 2 + 20) {
			n(0, j) = n_i;
		}
	}
	for (unsigned int t = 0; t < time_count; ++t) {
		for (unsigned int j = 0; j < pos_count; ++j) {
			if (j < pos_count / 2)
				Ca(t, j) = n_i;
			else
				Cd(t, j) = n_i;
		}
	}

	std::ofstream file("electric_field.txt");
	if (!file.is_open()) {
		std::cerr << "Could not open file electric_field.txt" << std::endl;
		return 1;
	}
	for (unsigned int t = 0; t < time_count - 1; ++t) {
		// We set boundary conditions: zero charges at the ends of our space
		n(t, -1) = 0;
		n(t, n.get_pos_count()) = 0;

		std::vector<double> electric_field;
		iter(n, p, Cd, Ca, electric_field, t);

		if (t % (time_count / 2000) == 0) {
			for (double e : electric_field) {
				file << e << " ";
			}
			file << std::endl;
		}
		/*if (t % 400 == 0) {
			double integ1 = 0.0;
			double integ2 = 0.0;
			for (unsigned int j = 0; j < pos_count; ++j) {
				integ1 += p(t, j);
				integ2 += n(t, j);
			}

			std::cout << integ1 << std::endl;
			std::cout << integ2 << std::endl << std::endl;
		}*/
	}

	Quantity Cd{ time_step, time_count, pos_step, pos_count };
	Quantity Ca{ time_step, time_count, pos_step, pos_count };
	for (unsigned int t = 0; t < time_count; ++t) {
		for (unsigned int j = 0; j < pos_count; ++j) {
			if (j < pos_count / 2) {
				Cd(0, j) = 0.001;
			}

			if (j >= pos_count / 2) {
				Ca(0, j) = 0.001;
			}
		}
	}


	Quantity rho{ time_end / static_cast<double>(time_count), time_count, pos_end / static_cast<double>(pos_count), pos_count };

	rho += p;
	rho += Cd;
	rho -= n;
	rho -= Ca;
	rho.print("data.dat", 2000);

	return 0;
}
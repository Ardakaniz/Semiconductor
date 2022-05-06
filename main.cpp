#include <iostream>
#include <vector>
#include <fstream>

/*
	get_pos_count() = X_max + 1
	get_time_count() = T_max + 1

	q(t, -1) = left boundary
	q(t, get_pos_count()) = right boundary
*/
class Quantity {
	public:
		Quantity(double time_step, unsigned int time_count, double pos_step, unsigned int pos_count, bool boundary_conditions = false) :
			m_time_step{ time_step },
			m_time_count{ time_count },
			m_pos_step{ pos_step },
			m_pos_count{ pos_count + 2 },
			m_boundary_conditions{ boundary_conditions }
		{
			for (unsigned int t = 0; t < m_time_count; ++t) {
				for (unsigned int p = 0; p < m_pos_count; ++p) {
					m_data.emplace_back(0.0);
				}
			}
		}
		Quantity(const Quantity& other) = default;
		Quantity(Quantity&& other) = default;
		Quantity& operator=(const Quantity& other) = default;
		Quantity& operator=(Quantity&& other) = default;

		Quantity& operator+=(const Quantity& other) {
			for (unsigned int t = 0; t < m_time_count; ++t) {
				for (unsigned int p = 0; p < m_pos_count; ++p) {
					m_data[t * m_pos_count + p] += other.m_data[t * m_pos_count + p];
				}
			}
			return *this;
		}

		Quantity& operator-=(const Quantity& other) {
			for (unsigned int t = 0; t < m_time_count; ++t) {
				for (unsigned int p = 0; p < m_pos_count; ++p) {
					m_data[t * m_pos_count + p] -= other.m_data[t * m_pos_count + p];
				}
			}
			return *this;
		}

		Quantity& operator*(double scalar) {
			for (unsigned int t = 0; t < m_time_count; ++t) {
				for (unsigned int p = 0; p < m_pos_count; ++p) {
					m_data[t * m_pos_count + p] *= scalar;
				}
			}
			return *this;
		}

		double get_time_step() const { return m_time_step; }
		double get_pos_step() const { return m_pos_step; }
		unsigned int get_time_count() const { return m_time_count; }
		unsigned int get_pos_count() const { return m_pos_count - 2; }

		double& operator()(int time, int pos) {
			int index = static_cast<int>(time * m_pos_count) + (pos + 1); // Because pos=-1 is associated with index 0

			if (!m_boundary_conditions) {
				// If we are on one of the edges and we dont have boundary conditions, we instead use the last modifiable index
				if (pos == -1)
					index++;
				else if (pos == get_pos_count())
					index--;
			}

			return m_data[static_cast<unsigned int>(index)];
		}

		double operator()(int time, int pos) const {
			int index = static_cast<int>(time * m_pos_count) + (pos + 1); // Because pos=-1 is associated with index 0

			if (!m_boundary_conditions) {
				// If we are on one of the edges and we dont have boundary conditions, we instead use the last modifiable index
				if (pos == -1)
					index++;
				else if (pos == get_pos_count())
					index--;
			}

			return m_data[static_cast<unsigned int>(index)];
		}

		void print(std::ostream& stream, unsigned int number_of_frame = 0, bool print_time = false) const {
			const unsigned int frame_freq = static_cast<unsigned int>(m_time_count / static_cast<double>(number_of_frame));
			for (unsigned int t = 0; t < m_time_count; ++t) {
				if (number_of_frame != 0 && t % frame_freq != 0)
					continue;

				if (print_time)
					stream << "[t = " << t * m_time_step << "]\t";

				for (unsigned int p = 1; p < m_pos_count - 1; ++p) {
					stream << m_data[t * m_pos_count + p] << " ";
				}
				stream << std::endl;
			}
		}

		void print(unsigned int number_of_frame = 0) const {
			print(std::cout, number_of_frame);
		}

		void print(const std::string& filename, unsigned int number_of_frame = 0) const {
			std::ofstream file(filename);
			if (!file.is_open()) {
				std::cerr << "Could not open file " << filename << std::endl;
				return;
			}

			print(file, number_of_frame);
		}

	private:
		const double m_time_step;
		const double m_pos_step;
		const unsigned int m_time_count;
		const unsigned int m_pos_count;
		const bool m_boundary_conditions;

		std::vector<double> m_data;
};

void iter(Quantity& n, Quantity& p, std::vector<double>& electric_field, unsigned int t) {
	const double K = 1.0; //* == \mu * e / \epsilon
	const double Cd = 0.0001;
	const double Ca = 0.0001;

	double qty_integral = n(t, 0) - p(t, 0);
	for (int j = 0; j < n.get_pos_count(); ++j) {
		const double dn = 0.5 * (n(t, j + 1) - n(t, j - 1));
		const double dp = 0.5 * (p(t, j + 1) - p(t, j - 1));

		double atoms_field = 0.0;
		if (j < n.get_pos_count() / 2)
			atoms_field = Cd * j;
		else
			atoms_field = Cd * ((n.get_pos_count() / 2) - 1) - Ca * (j - n.get_pos_count() / 2);

		const double E = (0.5 * (p(t, 0) + p(t, j) - n(t, 0) - n(t, j))) + qty_integral + atoms_field;
		
		electric_field.push_back(E * n.get_pos_step());
		n(t + 1, j) = n(t, j) + n.get_time_step() * K * (dn * E + (p(t, j) - n(t, j)) * n(t, j));
		p(t + 1, j) = p(t, j) - p.get_time_step() * K * (dp * E + (p(t, j) - n(t, j)) * p(t, j));

		qty_integral += p(t, j) - n(t, j);
	}
}

int main() {
	const double time_end = 5000;
	const unsigned int time_count = 50000;
	const double pos_end = 1.0;
	const unsigned int pos_count = 500;


	Quantity n{ time_end / static_cast<double>(time_count), time_count, pos_end / static_cast<double>(pos_count), pos_count };
	Quantity p{ time_end / static_cast<double>(time_count), time_count, pos_end / static_cast<double>(pos_count), pos_count };
	Quantity rho{ time_end / static_cast<double>(time_count), time_count, pos_end / static_cast<double>(pos_count), pos_count };

	for (unsigned int j = 0; j < pos_count; ++j) {
		n(0, j) = 0;
		p(0, j) = 0;

		if (j < pos_count / 2) {
			n(0, j) = 0.0001;
		}

		if (j >= pos_count / 2) {
			p(0, j) = 0.0001;
		}
	}

	std::ofstream file("electric_field.txt");
	if (!file.is_open()) {
		std::cerr << "Could not open file electric_field.txt" << std::endl;
		return 1;
	}
	for (unsigned int t = 0; t < time_count - 1; ++t) {
		std::vector<double> electric_field;
		iter(n, p, electric_field, t);

		if (t % 50 == 0) {
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

	rho += p;
	rho -= n;
	rho.print("data.dat", 1000);

	return 0;
}
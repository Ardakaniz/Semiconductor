#include <vector>
#include <string>
#include <ostream>

/*
	get_pos_count() = X_max + 1
	get_time_count() = T_max + 1

	q(t, -1) = left boundary
	q(t, get_pos_count()) = right boundary
*/
class Quantity {
	public:
		Quantity(double time_step, unsigned int time_count, double pos_step, unsigned int pos_count, bool boundary_conditions = false, bool fixed_time_count = false);
		Quantity(const Quantity& other) = default;
		Quantity(Quantity&& other) = default;
		Quantity& operator=(const Quantity& other) = default;
		Quantity& operator=(Quantity&& other) = default;

		Quantity& operator+=(const Quantity& other);
		Quantity& operator-=(const Quantity& other);
		Quantity& operator*(double scalar);

		double get_time_step() const { return m_time_step; }
		double get_pos_step() const { return m_pos_step; }
		unsigned int get_time_count() const { return m_time_count; }
		unsigned int get_pos_count() const { return m_pos_count - 2; }

		double& operator()(unsigned int time, int pos);
		double operator()(unsigned int time, int pos) const;

		void print(std::ostream& stream, unsigned int number_of_frame = 0, bool print_time = false) const;
		void print(unsigned int number_of_frame = 0) const;
		void print(const std::string& filename, unsigned int number_of_frame = 0) const;

	private:
		void add_time();
		void add_time(unsigned int count);

		const double m_time_step;
		const double m_pos_step;
		unsigned int m_time_count;
		const unsigned int m_pos_count;
		const bool m_boundary_conditions;
		const bool m_fixed_time_count;

		std::vector<double> m_data;
};
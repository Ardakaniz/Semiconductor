#include "Quantity.hpp"

#include <iostream>
#include <fstream>

Quantity::Quantity(double time_step, unsigned int time_count, double pos_step, unsigned int pos_count, bool boundary_conditions, bool fixed_time_count) :
	m_time_step{ time_step },
	m_time_count{ time_count },
	m_pos_step{ pos_step },
	m_pos_count{ pos_count + 2 },
	m_boundary_conditions{ boundary_conditions },
	m_fixed_time_count{ fixed_time_count }
{
	for (unsigned int t = 0; t < m_time_count; ++t) {
		for (unsigned int p = 0; p < m_pos_count; ++p) {
			m_data.emplace_back(0.0);
		}
	}
}

Quantity& Quantity::operator+=(const Quantity& other) {
	if (m_pos_count != other.m_pos_count)
		throw std::runtime_error("Quantity::operator+=: incompatible position count");
	if (m_time_count < other.m_time_count) {
		if (m_fixed_time_count)
			throw std::runtime_error("Quantity::operator+=: incompatible time count");
		else
			add_time(other.m_time_count - m_time_count);
	}

	for (unsigned int t = 0; t < other.m_time_count; ++t) {
		for (unsigned int p = 0; p < other.m_pos_count; ++p) {
			m_data[t * m_pos_count + p] += other.m_data[t * m_pos_count + p];
		}
	}
	return *this;
}

Quantity& Quantity::operator-=(const Quantity& other) {
	if (m_pos_count != other.m_pos_count)
		throw std::runtime_error("Quantity::operator+=: incompatible position count");
	if (m_time_count < other.m_time_count) {
		if (m_fixed_time_count)
			throw std::runtime_error("Quantity::operator+=: incompatible time count");
		else
			add_time(other.m_time_count - m_time_count);
	}

	for (unsigned int t = 0; t < other.m_time_count; ++t) {
		for (unsigned int p = 0; p < other.m_pos_count; ++p) {
			m_data[t * m_pos_count + p] -= other.m_data[t * m_pos_count + p];
		}
	}
	return *this;
}

Quantity& Quantity::operator*(double scalar) {
	for (unsigned int t = 0; t < m_time_count; ++t) {
		for (unsigned int p = 0; p < m_pos_count; ++p) {
			m_data[t * m_pos_count + p] *= scalar;
		}
	}
	return *this;
}

double& Quantity::operator()(unsigned int time, int pos) {
	if (!m_fixed_time_count && time >= m_time_count)
		add_time(m_time_count - time);

	if (time > m_time_count)
		throw std::runtime_error{"time index out of range"};
	if (pos < -1 || pos > static_cast<int>(get_pos_count()))
		throw std::runtime_error{"pos index out of range"};

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

double Quantity::operator()(unsigned int time, int pos) const {
	if (time > m_time_count)
		throw std::runtime_error{"time index out of range"};
	if (pos < -1 || pos > static_cast<int>(get_pos_count()))
		throw std::runtime_error{"pos index out of range"};
		
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

void Quantity::print(std::ostream& stream, unsigned int number_of_frame, bool print_time) const {
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

void Quantity::print(unsigned int number_of_frame) const {
	print(std::cout, number_of_frame);
}

void Quantity::print(const std::string& filename, unsigned int number_of_frame) const {
	std::ofstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Could not open file " << filename << std::endl;
		return;
	}

	print(file, number_of_frame);
}

void Quantity::add_time() {
	m_time_count++;
	for (unsigned int p = 0; p < m_pos_count; ++p) {
		m_data.emplace_back(0.0);
	}
}

void Quantity::add_time(unsigned int count) {
	for (unsigned int i = 0; i < count; ++i)
		add_time();
}
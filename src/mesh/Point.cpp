#include "Point.h"

Point::Point(const std::string &name,
			 const vecDouble &coordinates,
			 const double &lcar,
			 const bool &discretization)
{
	name_ = name;
	lcar_ = lcar;
	discretization_ = discretization;
	coordinates_.resize(3, 0.0);
	coordinates_[0] = coordinates[0];
	coordinates_[1] = coordinates[1];
	if (coordinates.size() == 2)
	{
		coordinates_[2] = 0.0;
	}
	else
	{
		coordinates_[2] = coordinates[2];
	}
}
Point::~Point() {}

std::string Point::getName()
{
	return name_;
}

std::string Point::getGmshCode()
{
	std::stringstream text;
	text << std::fixed;
	if (discretization_)
	{
		text << name_ << " = newp; Point(" << name_ << ") = {" << coordinates_[0] << ", " << coordinates_[1] << ", 0.0, " << lcar_
			 << "}; Physical Point('" << name_ << "') = {" << name_ << "};\n//\n";
	}
	else
	{
		text << name_ << " = newp; Point(" << name_ << ") = {" << coordinates_[0] << ", " << coordinates_[1] << ", 0.0, " << lcar_ << "};\n//\n";
	}
	return text.str();
}

void Point::addNodeToPoint(Node *node)
{
	pointNode_ = node;
}

Node *Point::getPointNode()
{
	return pointNode_;
}

void Point::getCoordinates(vecDouble &coord)
{
	coord[0] = coordinates_[0];
	coord[1] = coordinates_[1];
	coord[2] = coordinates_[2];
}

double Point::getlcar()
{
	return lcar_;
}

void Point::setlcar(const double &newlcar)
{
	lcar_ = newlcar;
}

void Point::setCoordinates(const vecDouble &newCoordinates)
{
	coordinates_[0] = newCoordinates[0];
	coordinates_[1] = newCoordinates[1];
	coordinates_[2] = newCoordinates[2];
}

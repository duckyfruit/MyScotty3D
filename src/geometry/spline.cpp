
#include "../geometry/spline.h"

template<typename T> T Spline<T>::at(float time) const {

	if(!any()) //if there are no knots, no splines!
	return T();
	
	std::set<float> pos = keys(); //get the keys

	int posS = int(pos.size()); //set the size to be posS
	std::vector<float> times;

	for(std::set<float>::iterator it=pos.begin(); it!=pos.end(); ++it) //put everything in a vector
	{
		float num = *it;
		times.push_back(num);
	}



	if((posS == 1))//if there is only one knot, return it!
	return knots.at(times[0]);

	if((has(time))) //if time is already there, no need to spline!
	return knots.at(time);

	if(time < times[0])
	return knots.at(times[0]); //if time is before the first knot, return the first knot

	if(time > times[posS- 1])
	return knots.at(times[posS - 1]); //if time is after the last knot, return the last knot

	//after this, size should be >= 2 and inside the bounds of the knots, therefore, we can interpolate!!
	int t0 = 0;
	int t1 = 0;
	int t2 = 0;
	int t3 = 0;
	
	T k1 = knots.at(times[0]); //set initial value
	T k0;
	T k2;
	T k3;
	T m0;
	T m1;

	

	for(int i = 0; i < posS; i++)
	{
		if(time > times[i]) //lower bound
		{
			t1 = i;
			k1 = knots.at(times[i]);
			

		}
	}

	t2 = t1 + 1; //upperbound 
	k2 = knots.at(times[t2]);

	if(t2 == (posS - 1)) //virtual knots
	{
		k3 = k2 + (k2 - k1);
		t3 = -1;
	}
	else //real knots
	{
		t3 = t2 + 1;
		k3 = knots.at(times[t3]);
	}

	if(t1 == 0) //virtual knots
	{
		k0 = k1 - (k2 - k1);
		t0 = -1;
	}
	else //real knots
	{
		t0 = t1 - 1;
		k0 = knots.at(times[t0]);
	}

	//printf("%d", t0);

	float interval = times[t2] - times[t1]; //interval range of time!
	//need to divide by the time range 
	if(t0 == -1)
	m0 = ( k2 - k0 ) /  ( times[t2]);
	else
	m0 = ( k2 - k0 ) /  ( times[t2] - times[t0] );

	if(t3 == -1)
	m1 = ( k3 - k1 ) /  (-times[t1] ); 
	else
	m1 = ( k3 - k1 ) /  ( times[t3] - times[t1] ); 


	


	//go throught the keys list and find out which knots time is between
	//check if there is a second knot before this knot


	

	// A4T1b: Evaluate a Catumull-Rom spline

	// Given a time, find the nearest positions & tangent values
	// defined by the control point map.

	// Transform them for use with cubic_unit_spline

	// Be wary of edge cases! What if time is before the first knot,
	// before the second knot, etc...

	return cubic_unit_spline(time/interval, (k1 - k1), (k2 - k1), m0, m1);
}

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1) {
	
	float timeS = pow(time, 2.0f);
	float timeC = pow(time, 3.0f);

	float h00 =  2 * timeC + 3* timeS + 1;
	float h01 = timeC - 2* timeS + time;
	float h10 = -2 * timeC + 3 * timeS;
	float h11 = timeC - timeS;

	T interp = h00 * position0 + h10 * tangent0 + h01 * position1 + h11 * tangent1; 

	// A4T1a: Hermite Curve over the unit interval

	// Given time in [0,1] compute the cubic spline coefficients and use them to compute
	// the interpolated value at time 'time' based on the positions & tangents

	// Note that Spline is parameterized on type T, which allows us to create splines over
	// any type that supports the * and + operators.

	return  interp;
}

template class Spline<float>;
template class Spline<double>;
template class Spline<Vec4>;
template class Spline<Vec3>;
template class Spline<Vec2>;
template class Spline<Mat4>;
template class Spline<Spectrum>;

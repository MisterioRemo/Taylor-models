#include "stdafx.h"
#include "CppUnitTest.h"
#include "../TaylorModel/interval.h"
#include <algorithm>
using namespace Microsoft::VisualStudio::CppUnitTestFramework;


namespace intervalClassTest
{		

	TEST_CLASS(UnitTest1)
	{
	public:
		
		TEST_METHOD(TestMethodConstructor)
		{
			interval<double> i1(1.25, 0.6);
			interval<double> i2(i1);
			Assert::IsTrue(i1 == i2);
		}

		TEST_METHOD(TestMethodPlus1)
		{
			interval<double> i1(1.2, 4.6);
			interval<double> i2(3.8, 2.6);

			i1 += i2;
			Assert::IsTrue(i1 == interval<double>(1.2 + 3.8, 4.6 + 2.6));
		}

		TEST_METHOD(TestMethodPlus2)
		{
			interval<double> i1(1.2, 4.6);
			interval<double> i2(3.8, 2.6);

			interval<double> i3 = i1 + i2;
			Assert::IsTrue(i3 == interval<double>(1.2 + 3.8, 4.6 + 2.6));
		}

		TEST_METHOD(TestMethodMinus1)
		{
			interval<double> i1(1.2, 4.6);
			interval<double> i2(3.8, 2.6);

			i1 -= i2;
			Assert::IsTrue(i1 == interval<double>(1.2 - 2.6, 4.6 - 3.8));
		}

		TEST_METHOD(TestMethodMinus2)
		{
			interval<double> i1(1.2, 4.6);
			interval<double> i2(3.8, 2.6);

			interval<double> i3 = i1 - i2;
			Assert::IsTrue(i3 == interval<double>(1.2 - 2.6, 4.6 - 3.8));
		}

		TEST_METHOD(TestMethodMult1)
		{
			interval<double> i1(1.2, 4.6);
			interval<double> i2(3.8, 2.6);
			i1 *= i2;

			double arr[] = { 1.2*3.8, 4.6*2.6, 1.2*2.6, 4.6*3.8 };
			double begin = *std::min_element(arr, arr + 4);
			double end = *std::max_element(arr, arr + 4);
			Assert::IsTrue(i1 == interval<double>(begin, end));
		}

		TEST_METHOD(TestMethodMult2)
		{
			interval<double> i1(1.2, 4.6);
			interval<double> i2(3.8, 2.6);
			interval<double> i3 = i1 * i2;

			double arr[] = { 1.2*3.8, 4.6*2.6, 1.2*2.6, 4.6*3.8 };
			double begin = *std::min_element(arr, arr + 4);
			double end = *std::max_element(arr, arr + 4);
			Assert::IsTrue(i3 == interval<double>(begin, end));
		}

		TEST_METHOD(TestMethodMult3)
		{
			interval<double> i1(1.2, 4.6);
			double t = 3.0045;
			i1 *= t;

			Assert::IsTrue(i1 == interval<double>(1.2 * t, 4.6 * t));
		}

		TEST_METHOD(TestMethodMult4)
		{
			double t = 3.0045;
			interval<double> i1(1.2, 4.6);
			interval<double> i2 = i1 * t;

			Assert::IsTrue(i2 == interval<double>(1.2 * t, 4.6 * t));
		}

		TEST_METHOD(TestMethodDiv1)
		{
			bool exceptionThrown = false;
			interval<double> i1(1, 2);
			interval<double> i2(1, 0);
			try {
				interval<double> i3 = i1 / i2;
			}
			catch (interval<double>::divideByZero &ex) {
				exceptionThrown = true;
			}

			Assert::IsTrue(exceptionThrown);
		}

	};
}
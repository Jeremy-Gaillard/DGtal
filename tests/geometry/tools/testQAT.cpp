/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file testPreimage.cpp
 * @ingroup Tests
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 *
 * @date 2010/07/02
 *
 * This file is part of the DGtal library
 */

/**
 * @brief Aim: simple test of Preimage2D
 * @see testGeometricalDSS.cpp
 */

#include <fstream>
#include <iostream>



#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/boards/Board2D.h"

#include "DGtal/base/Common.h"
#include "DGtal/kernel/SimpleMatrix.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/geometry/tools/NaiveQAT.h"
#include "DGtal/geometry/tools/NearestNeighborQAT.h"

#include "DGtal/images/CConstImage.h"

//#include "DGtal/io/readers/MagickReader.h"



#include "ConfigTest.h"


using namespace DGtal;
using namespace Z2i;

typedef ImageContainerBySTLVector<Z2i::Domain, int> Image;


void testNearestNeighborQAT(const Image & image, const SimpleMatrix<int, 2, 2> & Mat, const int & omega, const Point & vect)
{
  //BOOST_CONCEPT_ASSERT(( CConstImage<NearestNeighborQAT<Image> > ));
  
  NearestNeighborQAT<Image> QAT( Mat, omega, vect, 0 );
  QAT.setImage(image);
  Domain newDomain = QAT.domain();
  
  Board2D board;
  board << SetMode( newDomain.className(), "Paving" )
    << newDomain
    << SetMode( newDomain.lowerBound().className(), "Paving" );
  string specificStyle = newDomain.lowerBound().className() + "/Paving";
  
  GradientColorMap<int> cmap_grad( 1, 20 );
  cmap_grad.addColor( Color( 50, 50, 255 ) );
  cmap_grad.addColor( Color( 255, 0, 0 ) );
  cmap_grad.addColor( Color( 255, 255, 10 ) );
  
  typename NearestNeighborQAT<Image>::ConstRange::ConstIterator ite = QAT.constRange().end();
  for ( typename NearestNeighborQAT<Image>::ConstRange::ConstIterator it = QAT.constRange().begin();
       it != ite; it++ )
  {
    cout << *it << " ";
  }
  cout << endl;
  
  for ( typename Domain::Iterator it = newDomain.begin(); it != newDomain.end(); it++ )	
  {
    int value = QAT(*it);
    if( value == 0)
      board << CustomStyle( specificStyle,
	new CustomColors( Color::Black,
	Color::Black ) )
	<< *it;
    else
      board << CustomStyle( specificStyle,
	  new CustomColors( Color::Black,
	  cmap_grad( value ) ) )
	  << *it;
  }
  
  board.saveEPS("testNNQAT.eps");
}

void testNaiveQAT(const Image & image, const SimpleMatrix<int, 2, 2> & Mat, const int & omega, const Point & vect)
{
  //BOOST_CONCEPT_ASSERT(( CConstImage<NaiveQAT<Image> > ));
  
  NaiveQAT<Image> QAT( Mat, omega, vect, 0 );
  QAT.setImage(image);
  Domain newDomain = QAT.domain();
  
  Board2D board;
  board << SetMode( newDomain.className(), "Paving" )
    << newDomain
    << SetMode( newDomain.lowerBound().className(), "Paving" );
  string specificStyle = newDomain.lowerBound().className() + "/Paving";
  
  GradientColorMap<int> cmap_grad( 1, 20 );
  cmap_grad.addColor( Color( 50, 50, 255 ) );
  cmap_grad.addColor( Color( 255, 0, 0 ) );
  cmap_grad.addColor( Color( 255, 255, 10 ) );
  
  typename NaiveQAT<Image>::ConstRange::ConstIterator ite = QAT.constRange().end();
  for ( typename NaiveQAT<Image>::ConstRange::ConstIterator it = QAT.constRange().begin();
       it != ite; it++ )
  {
    cout << *it << " ";
  }
  cout << endl;
  
  for ( typename Domain::Iterator it = newDomain.begin(); it != newDomain.end(); it++ )	
  {
    int value = QAT(*it);
    if( value == 0)
      board << CustomStyle( specificStyle,
	new CustomColors( Color::Black,
	Color::Black ) )
	<< *it;
    else
      board << CustomStyle( specificStyle,
	  new CustomColors( Color::Black,
	  cmap_grad( value ) ) )
	  << *it;
  }
  
  board.saveEPS("testNaiveQAT.eps");
}


int main()
{
  Point p1(0,0);
  Point p2(4,4);
 
  Domain domain(p1, p2);
  
  
  Image image(domain);
  image.setValue(Point(0,0), 1); image.setValue(Point(1,0), 7); image.setValue(Point(2,0), 8); image.setValue(Point(3,0), 9); image.setValue(Point(4,0), 1);
  image.setValue(Point(0,1), 13); image.setValue(Point(1,1), 1); image.setValue(Point(2,1), 2); image.setValue(Point(3,1), 1); image.setValue(Point(4,1), 11);
  image.setValue(Point(0,2), 12); image.setValue(Point(1,2), 1); image.setValue(Point(2,2), 2); image.setValue(Point(3,2), 1); image.setValue(Point(4,2), 12);
  image.setValue(Point(0,3), 11); image.setValue(Point(1,3), 1); image.setValue(Point(2,3), 3); image.setValue(Point(3,3), 2); image.setValue(Point(4,3), 13);
  image.setValue(Point(0,4), 1); image.setValue(Point(1,4), 9); image.setValue(Point(2,4), 8); image.setValue(Point(3,4), 7); image.setValue(Point(4,4), 1);

  //Image image = MagickReader<Image>::importImage("inQAT.png");
  //Domain domain = image.domain();
  

  SimpleMatrix<int, 2, 2> Mat;
  Mat.setComponent(0, 0, 5);
  Mat.setComponent(0, 1, -2);
  Mat.setComponent(1, 0, 1);
  Mat.setComponent(1, 1, 4);
//   Mat.setComponent(0, 0, 4);
//   Mat.setComponent(0, 1, 1);
//   Mat.setComponent(1, 0, 1);
//   Mat.setComponent(1, 1, 1);
  int omega = 2;
  Point vect(0, 0);
  
  testNaiveQAT(image, Mat, omega, vect);
  
  testNearestNeighborQAT(image, Mat, omega, vect);
  
}

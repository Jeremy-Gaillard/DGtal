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

#pragma once

/**
 * @file NaiveQAT.h
 * @author Jérémy Gaillard (\c jeremy.gaillard@insa-lyon.fr )
 * Institut National des Sciences Appliquées - INSA, France
 *
 * @date 2012/07/10
 *
 * @brief Header file for module NaiveQAT.h
 *
 * This file is part of the DGtal library.
 */

#if defined(NaiveQAT_RECURSES)
#error Recursive header files inclusion detected in Preimage2D.h
#else // defined(NaiveQAT_RECURSES)
/** Prevents recursive inclusion of headers. */
#define NaiveQAT_RECURSES

#if !defined NaiveQAT_h
/** Prevents repeated inclusion of headers. */
#define NaiveQAT_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/SimpleMatrix.h"
#include "DGtal/kernel/sets/DigitalSetSelector.h"
#include "DGtal/base/CowPtr.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class NaiveQAT
  /**
   * @brief Aim: Computes the preimage of the 2D Euclidean shapes 
   * crossing a sequence of n straigth segments in O(n),
   * with the algorithm of O'Rourke (1981). 
   *
   * \b Model of CConstImage
   *
   *
   * @tparam TImageContainer any model of CImage
   */
  template <typename TImageContainer>
  class NaiveQAT
  {
    //BOOST_CONCEPT_ASSERT(( CImage<TImageContainer> ));


    // ----------------------- Types ------------------------------
  public:
    typedef TImageContainer ImageContainer;
    typedef typename TImageContainer::Domain Domain;
    typedef typename TImageContainer::Point Point;
    typedef typename TImageContainer::Value Value;
    typedef typename TImageContainer::ConstRange ConstRange;
    typedef CowPtr<TImageContainer> ImagePointer;
    
    typedef typename Domain::Dimension Dimension;
    
    static const Dimension dimension = Domain::dimension;
    
    typedef SimpleMatrix<Value, dimension, dimension> Matrix;
    typedef typename Matrix::ColumnVector Vector;	// Same as point
    typedef typename DigitalSetSelector
	< Domain, SMALL_DS+HIGH_ITER_DS >
	::Type Paving;
    
    

  private:

    
    


    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     */
    NaiveQAT( const Matrix & M, const Value & omega, const Vector & V );
    
    /**
     * Destructor. Does nothing.
     */
    ~NaiveQAT();

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    NaiveQAT( const NaiveQAT & other );
    
    
    /**
     * Returns a reference to the transformed image domain.
     *
     * @return a reference to the domain.
     */
    const Domain & domain() const
    {
      return myImage->domain();
    }
    
    /**
     * Returns the range of the underlying image
     * to iterate over its values
     *
     * @return a range.
     */
    ConstRange constRange() const;
    
    /**
     * Get the value of the transformed image at a given 
     * position given by a point.
     *
     * @pre the point must be in the transformed image domain
     *
     * @param aPoint the point.
     * @return the value at aPoint.
     */
    Value operator()( const Point & aPoint ) const
    {
      return (*myImage)(aPoint);
    }
    
    
    /**
     * Applies the QAT to an image. The result is stored in @a myImage.
     * 
     * @param image the image to be transformed
     */
    void transformImage( const ImageContainer & image );
    

    

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;


    // ------------------------- Protected Datas ------------------------------
  protected:
    Matrix myM;
    Vector myV;
    Value myOmega;
    /**
     * Inverted components
     */
    Matrix myMInv;
    Vector myVInv;
    Value myOmegaInv;
    
    ImagePointer myImage;
    
    std::vector<Paving> myPavings;
    std::vector<Paving> myPavingsRemainder;
    
    
    // ------------------------- Private Datas --------------------------------
  private:
    

    // ------------------------- Hidden services ------------------------------
  protected:
    NaiveQAT();
    
    /**
      * computes the image of a point
      *
      * @param p the initial point
      *
      * @return the image of the point
      */
    const Point calculate ( const Point & p ) const;
    
    /**
      * computes the minimum and maximum coordinates in the final image : min_x, min_y, max_x, max_y
      *
      * @param domain the initial image's domain
      *
      * @return the matrix : (min_i, min_j, max_i, max_j)
      */
    Domain getImageBound ( const Domain domain );
    
    Value backwardColorLinear( const Point & Pp, const ImageContainer & image);
    
    void recursiveColor( double & val, double & div, const std::vector<double> & x,
			 const std::vector<double> & r, const std::vector<double> & l,
			 const ImageContainer & image, int start, const std::vector<bool> & sequence ) const;
    
    void getAllVertices( std::vector<Point> & vertices, const Dimension & start, 
		const Point & startPoint, const Point & endPoint ) const;
		
    void recursiveTransform( const Point & lowerBound, const Point & upperBound,
			     Dimension dim, std::vector<typename Point::Component> & components,
			     const std::vector<Point> & incr, Point & Pp, const ImageContainer & image );

    void inverse();

  private:

    

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class NaiveQAT


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/geometry/tools/NaiveQAT.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined NaiveQAT_h

#undef NaiveQAT_RECURSES
#endif // else defined(NaiveQAT_RECURSES)

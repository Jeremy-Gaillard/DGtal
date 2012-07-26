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
 * @file FastQAT.h
 * @author Jérémy Gaillard (\c jeremy.gaillard@insa-lyon.fr )
 * Institut National des Sciences Appliquées - INSA, France
 *
 * @date 2012/07/10
 *
 * @brief Header file for module FastQAT.h
 *
 * This file is part of the DGtal library.
 */

#if defined(FastQAT_RECURSES)
#error Recursive header files inclusion detected in Preimage2D.h
#else // defined(FastQAT_RECURSES)
/** Prevents recursive inclusion of headers. */
#define FastQAT_RECURSES

#if !defined FastQAT_h
/** Prevents repeated inclusion of headers. */
#define FastQAT_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/sets/DigitalSetSelector.h"
#include "DGtal/kernel/SimpleMatrix.h"
#include "DGtal/base/CowPtr.h"
#include "DGtal/images/DefaultConstImageRange.h"
#include "DGtal/images/CConstImage.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class FastQAT
  /**
   * @brief Aim : computes a quasi-affine transformation (M+V)/omega 
   * of an image using a fast algorithm
   *
   * \b Model of CConstImage
   *
   *
   * @tparam TImageContainer any model of CConstImage
   */
  template<typename TImageContainer, Dimension d = TImageContainer::Domain::dimension>
  class FastQAT
  {
    FastQAT()
    {
      assert( false );
    }
  };
  
  
  template <typename TImageContainer>
  class FastQAT<TImageContainer, 2>
  {
    BOOST_CONCEPT_ASSERT(( CConstImage<TImageContainer> ));


    // ----------------------- Types ------------------------------
  public:
    typedef FastQAT< TImageContainer > Self;
    typedef TImageContainer ImageContainer;
    typedef typename TImageContainer::Domain Domain;
    typedef typename TImageContainer::Point Point;
    typedef typename TImageContainer::Value Value;
    typedef CowPtr<const TImageContainer> ImagePointer;
    
    typedef DefaultConstImageRange< Self > ConstRange;
    
    typedef typename DigitalSetSelector
	< Domain, SMALL_DS+HIGH_ITER_DS >
	::Type Paving;

    typedef typename Domain::Dimension Dimension;
    
    static const Dimension dimension = Domain::dimension;

    typedef SimpleMatrix<Value, dimension, dimension> Matrix;
    typedef typename Matrix::ColumnVector Vector;	// Same as point
    
    
  private:

    
    


    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * 
     * @param M transformation matrix
     * 
     * @param omega
     * 
     * @param V translation vector
     * 
     * @param defaultValue the value affected to the points outside the 
     * transformed image
     */
    FastQAT( const Matrix & M, const Value & omega, const Vector & V, const Value & defaultValue, const ImageContainer & image );

    /**
     * Destructor. Does nothing.
     */
    ~FastQAT() { }

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    FastQAT( const FastQAT & other )
    { 
      myM = other.myM;
      myV = other.myV;
      myDefaultValue = other.myDefaultValue;
      myImage = other.myImage;
      myOmega = other.myOmega;
      myDomain = other.myDomain;
      isContracting = other.isContracting;
      pavings = other.pavings;
      pavingsRemainder = other.pavingsRemainder;
      myMInv = other.myMInv;
      myVInv = other.myVInv;
      myOmegaInv = other.myOmegaInv;
      compute();	// TODO : replace by a0=other.a0 ...
    }
    
    /**
     * Chooses the image that will be transformed by the QAT.
     * 
     * @param image the image to be transformed
     */
    //void setImage( const ImageContainer & image );
    
    /**
     * Returns a reference to the transformed image domain.
     *
     * @return a reference to the domain.
     */
    const Domain & domain() const
    {
      return myDomain;
    }
    
    /**
     * Returns the range of the underlying image
     * to iterate over its values
     *
     * @return a range.
     */
    ConstRange constRange() const
    {
      return ConstRange( *this );
    }
    
    /**
     * Get the value of the transformed image at a given 
     * position given by a Point.
     *
     * @pre the point must be in the transformed image domain
     *
     * @param aPoint the point.
     * @return the value at aPoint.
     */
    Value operator()( const Point & aPoint ) const;

    

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const
    {
      myImage->selfDisplay();
    }


    // ------------------------- Protected Datas ------------------------------
  protected:
    
    Value myDefaultValue;
    /**
     * Transformed image domain
     */
    Domain myDomain;
    
    bool isContracting;
    
    /**
     * QAT parameters
     */
    Matrix myM;
    Vector myV;
    Value myOmega;
    /**
     * Inverted components
     */
    Matrix myMInv;
    Vector myVInv;
    Value myOmegaInv;
    
    /**
     * Transformed image
     */
    ImagePointer myImage;
    
    std::vector<Paving> pavings;
    std::vector<Paving> pavingsRemainder;
    
    Value a0, b0, c0, d0, e0, f0;
    Value a1, b1, c1, d1;
    Value b10, b11;
    Matrix H, H0;
    Value alpha0, alpha1, beta1;
    Vector U0, U1;
    
    // ------------------------- Private Datas --------------------------------
  private:
    

    // ------------------------- Hidden services ------------------------------
  protected:
    FastQAT();
    
    /**
      * Computes the domain of the final image
      *
      * @param domain the initial image's domain
      *
      * @return the final image's domain
      */
    Domain getImageBound ( const Domain domain ) const;
    
    /**
      * Computes the image of a point
      *
      * @param p the initial point
      *
      * @return the image of the point
      */
    const Point calculate ( const Point & p ) const;
    
    /**
      * Computes the antecedent of a point
      *
      * @param p the initial point
      *
      * @return the image of the point
      */
    const Point calculateInv ( const Point & p ) const;
    
    /**
      * computes the remainder of a point
      *
      * @param p the initial point
      *
      * @return the image of the point
      */
    const Point calculateRemainderInv ( const Point p ) const;

    void setPavingRemainder ( const Point I, const Point Rem, const Point P );
    
    void determinePavingsWithRemainders();
    
    Value pavingValue( const Paving & P, const Point & v ) const;
    
    /**
     * Computes the inverse of the QAT
     */
    void inverse();
    
    void compute();
    
    


  private:

    

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class FastQAT


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/geometry/tools/FastQAT.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined FastQAT_h

#undef FastQAT_RECURSES
#endif // else defined(FastQAT_RECURSES)

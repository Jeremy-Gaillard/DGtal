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
 * @file NearestNeighborQAT.h
 * @author Jérémy Gaillard (\c jeremy.gaillard@insa-lyon.fr )
 * Institut National des Sciences Appliquées - INSA, France
 *
 * @date 2012/07/10
 *
 * @brief Header file for module NearestNeighborQAT.h
 *
 * This file is part of the DGtal library.
 */

#if defined(NearestNeighborQAT_RECURSES)
#error Recursive header files inclusion detected in Preimage2D.h
#else // defined(NearestNeighborQAT_RECURSES)
/** Prevents recursive inclusion of headers. */
#define NearestNeighborQAT_RECURSES

#if !defined NearestNeighborQAT_h
/** Prevents repeated inclusion of headers. */
#define NearestNeighborQAT_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/sets/DigitalSetSelector.h"
#include "DGtal/kernel/SimpleMatrix.h"
#include "DGtal/base/CowPtr.h"
#include "DGtal/images/CConstImage.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class NearestNeighborQAT
  /**
   * @brief Aim : computes a quasi-affine transformation (M+V)/omega 
   * of an image using the nearest neighbor algorithm
   *
   * \b Model of CConstImage
   *
   *
   * @tparam TImageContainer any model of CConstImage
   */
  template <typename TImageContainer>
  class NearestNeighborQAT
  {
    BOOST_CONCEPT_ASSERT(( CConstImage<TImageContainer> ));


    // ----------------------- Types ------------------------------
  public:
    typedef TImageContainer ImageContainer;
    typedef typename TImageContainer::Domain Domain;
    typedef typename TImageContainer::Point Point;
    typedef typename TImageContainer::Value Value;
    typedef typename TImageContainer::ConstRange ConstRange;
    typedef CowPtr<TImageContainer> ImagePointer;
    
    typedef typename DigitalSetSelector
	< Domain, SMALL_DS+HIGH_ITER_DS >
	::Type Paving;
	
    typedef typename Domain::Dimension Dimension;
    
    static const Dimension dimension = Domain::dimension;
    

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
     */
    NearestNeighborQAT( const Matrix & M, const Value & omega, const Vector & V );

    /**
     * Destructor. Does nothing.
     */
    ~NearestNeighborQAT();

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    NearestNeighborQAT( const NearestNeighborQAT & other );
    
    
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
    ConstRange constRange() const
    {
      return myImage->constRange();
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
    
    // ------------------------- Private Datas --------------------------------
  private:
    

    // ------------------------- Hidden services ------------------------------
  protected:
    NearestNeighborQAT();
    
    /**
      * Computes the image of a point
      *
      * @param p the initial point
      *
      * @return the image of the point
      */
    const Point calculate ( const Point & p ) const;
    
    /**
     * Computes the inverse of the QAT
     */
    void inverse();


  private:

    

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class NearestNeighborQAT


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/geometry/tools/NearestNeighborQAT.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined NearestNeighborQAT_h

#undef NearestNeighborQAT_RECURSES
#endif // else defined(NearestNeighborQAT_RECURSES)

#include <fantom/algorithm.hpp>
#include <fantom/graphics.hpp>
#include <fantom/dataset.hpp>
#include <fantom/register.hpp>
 
// needed for BoundinSphere-Calculation and normal calculation
#include <fantom-plugins/utils/Graphics/HelperFunctions.hpp>
#include <fantom-plugins/utils/Graphics/ObjectRenderer.hpp>

#include <map>
#include <climits>

#ifndef M_PI
#define     M_PI        3.14159265358979323846
#endif

#define     EQN_EPS     1e-9
#define	    IsZero(x)	((x) > -EQN_EPS && (x) < EQN_EPS)

#ifdef NOCBRT
#define     cbrt(x)     ((x) > 0.0 ? pow((double)(x), 1.0/3.0) : \
                          ((x) < 0.0 ? -pow((double)-(x), 1.0/3.0) : 0.0))
#endif

 
using namespace fantom;
 
namespace
{
 
    class MarchingDiamondsAlgorithm : public VisAlgorithm
    {

    private:
        typedef std::pair< const size_t, const size_t > EdgeKey;

        struct Tetrahedron {
            std::array< size_t, 4 > vertices;
            uint8_t splitLevel;
        };

        // Store mutable copy of the input grid's vertices and cell information
        std::vector< Vector3 > gridPoints;
        std::vector< Tetrahedron > gridCells;

        // Each entry corresponds to an edge and contains a list of cell indices that are incident to this edge
        std::map< EdgeKey, std::vector< size_t > > incidentCellsOfEdge;
        // 
        std::map< EdgeKey, Vector3 > edgeIntersections;

        std::unique_ptr< FieldEvaluator< 3, Scalar > > evaluator;


        // Helper function to find the real-valued solutions of a cubic polynomial
        // source: https://github.com/erich666/GraphicsGems/blob/master/gems/Roots3And4.c
        int solveCubicPolynomial( double c[ 4 ], double s[ 3 ] )
        {
            int     i, num;
            double  sub;
            double  A, B, C;
            double  sq_A, p, q;
            double  cb_p, D;

            /* normal form: x^3 + Ax^2 + Bx + C = 0 */

            A = c[ 2 ] / c[ 3 ];
            B = c[ 1 ] / c[ 3 ];
            C = c[ 0 ] / c[ 3 ];

            /*  substitute x = y - A/3 to eliminate quadric term:
            x^3 +px + q = 0 */

            sq_A = A * A;
            p = 1.0/3 * (- 1.0/3 * sq_A + B);
            q = 1.0/2 * (2.0/27 * A * sq_A - 1.0/3 * A * B + C);

            /* use Cardano's formula */

            cb_p = p * p * p;
            D = q * q + cb_p;

            if ( IsZero(D) ) {
                if ( IsZero(q) ) {  // one triple solution
                    s[ 0 ] = 0;
                    num = 1;
                } else {  // one single and one double solution
                    double u = cbrt(-q);
                    s[ 0 ] = 2 * u;
                    s[ 1 ] = - u;
                    num = 2;
                }
            } else if (D < 0) {  // Casus irreducibilis: three real solutions
                double phi = 1.0/3 * acos(-q / sqrt(-cb_p));
                double t = 2 * sqrt(-p);

                s[ 0 ] =   t * cos(phi);
                s[ 1 ] = - t * cos(phi + M_PI / 3);
                s[ 2 ] = - t * cos(phi - M_PI / 3);
                num = 3;
            } else {  // one real solution
                double sqrt_D = sqrt(D);
                double u = cbrt(sqrt_D - q);
                double v = - cbrt(sqrt_D + q);

                s[ 0 ] = u + v;
                num = 1;
            }

            /* resubstitute */
            sub = 1.0 / 3 * A;
            for ( i = 0; i < num; ++i ) s[ i ] -= sub;

            return num;
        }

    public:
        struct Options : public VisAlgorithm::Options
        {
            Options( fantom::Options::Control& control )
                : VisAlgorithm::Options( control )
            {
                add< Field< 3, Scalar > >( "Field", "Scalar field defined on an unstructured 3-dimensional grid" );
                add< double >( "Iso value", "", 1.0 );
                add< size_t >( "Max splits", "", 2 );
                add< bool >( "No diamonds?", "", false );
            }
        };

        struct VisOutputs : public VisAlgorithm::VisOutputs
        {
            
            VisOutputs( fantom::VisOutputs::Control& control )
                : VisAlgorithm::VisOutputs( control )
            {
                addGraphics( "Iso surface" );
            }
        };

        MarchingDiamondsAlgorithm( InitData& data )
            : VisAlgorithm( data )
        {}

        EdgeKey getEdgeKey( size_t pointIndex1, size_t pointIndex2 )
        {
            return std::make_pair( std::min( pointIndex1, pointIndex2 ), std::max( pointIndex1, pointIndex2 ) );
        }


        bool checkEdge( size_t pointIndex1, size_t pointIndex2 )
        {
            if ( incidentCellsOfEdge.find( getEdgeKey ( pointIndex1, pointIndex2 ) ) != incidentCellsOfEdge.end() ) {
                return true;
            }
            return false;
        }

        double getFieldValue( Vector3 position )
        {
            // Get the interpolated field value at a given position
            if ( !evaluator->reset( position ) ) {
                // If the position is out of the domain, an exception is thrown
                throw std::domain_error( "Cannot interpolate because point is out of the domain!" );
            }
            return evaluator->value()[ 0 ];
        }


        bool getDiamondPoints( EdgeKey edgeKey, std::vector< size_t >& diamondPoints )
        {
            // Get set of points from the tetrahedra cells that are not the endpoints of the incident edge e
            // They will potentially form the ring of the diamond
            std::set< size_t > ringPoints;
            for ( size_t cellIndex : incidentCellsOfEdge[ edgeKey ] ) {
                Tetrahedron cell = gridCells[ cellIndex ];
                
                for ( size_t vertex : cell.vertices ) {
                    if ( vertex == edgeKey.first || vertex == edgeKey.second ) {
                        // Cell point is part of the incident edge and must be skipped
                        continue;
                    }
                    // Insert to set, values are unique so there will be no duplicates
                    ringPoints.insert( vertex );
                }
            }

            // Count of ring points d_i (i = 0, 1, ..., k - 1) must be equal to the number of tetrahedra k,
            // otherwise the tetrahedra do not form a diamond and the incident edge can be skipped 
            if ( ringPoints.size() != incidentCellsOfEdge[ edgeKey ].size() ) {
                return false;
            }

            const size_t K = ringPoints.size();

            // Initialize diamondPoints list with first entry of the ringPoints set 
            auto ringIterator = ringPoints.begin();
            diamondPoints = { *ringIterator };
            ringPoints.erase( ringIterator );

            // Add points of the diamond's ring to the list `diamondPoints` such that subsequent points are always
            // connected with an edge
            while ( diamondPoints.size() < K ) {
                bool foundPoint = false;
                for ( ringIterator = ringPoints.begin(); ringIterator != ringPoints.end(); ringIterator++ ) {
                    // Get last point of the list
                    size_t lastPointIndex = diamondPoints.back();

                    if ( checkEdge( lastPointIndex, *ringIterator ) ) {
                        // If edge exists we add the current ring point to the diamond point list and remove it from
                        // the set `ringPoints`
                        diamondPoints.push_back( *ringIterator );
                        ringPoints.erase( ringIterator );
                        // Set flag to indicate that a suitable point was found
                        foundPoint = true;
                        break;
                    }
                }

                // If no suitable point was found, the ring is incomplete and exit the loop early
                if ( !foundPoint ) break;
            }

            // At last, we need to check whether an edge exists between the first and the last point
            // This maybe not necessary because of previous checks, but let's just be safe here
            if ( diamondPoints.size() < K || !checkEdge( diamondPoints.front(), diamondPoints.back() ) ) {
                return false;
            }

            // Add endpoints of incident edge to list of diamond points (d_k and d_k+1)
            diamondPoints.push_back( edgeKey.first );
            diamondPoints.push_back( edgeKey.second );
            // Shrink capacity of vector to current size to free unused memory
            diamondPoints.shrink_to_fit();
            return true;
        }

 
        virtual void execute( const Algorithm::Options& options, const volatile bool& /*abortFlag*/ ) override
        {
            std::shared_ptr< const Function< Scalar > > function = options.get< Function< Scalar > >( "Field" );
            std::shared_ptr< const Field< 3, Scalar > > field = options.get< Field< 3, Scalar > >( "Field" );

            // If there is no input, do nothing
            if( !field ) {
                debugLog() << "Input field is missing." << std::endl;
                return;
            }
 
            evaluator = field->makeEvaluator();

            std::shared_ptr< const Grid< 3 > > grid = std::dynamic_pointer_cast< const Grid< 3 > >( function->domain() );
            if ( !grid ) {
                throw std::logic_error( "Wrong type of grid!" );
            }

            const double iso = options.get< double >( "Iso value" );
            const size_t maxSplits = options.get< size_t >( "Max splits" );
            const bool noDiamonds = options.get< bool >( "No diamonds?" );
            

            //==========================================================================================================
            // STEP 1: Create auxiliary data structure that represent the relationship of edges and tetrahedra cells
            //==========================================================================================================

            // Reset maps between runs
            incidentCellsOfEdge.clear();
            edgeIntersections.clear();

            // Reset vector lists as well
            gridPoints.clear();
            gridCells.clear();


            // Make mutable copy of points array
            const ValueArray< Point3 >& points = grid->points();
            for ( size_t i = 0; i < grid->numPoints(); ++i ) {
                gridPoints.push_back( (Vector3)points[ i ] );
            }

            // Fill list `gridCells` with Tetrahedron cells from grid
            for ( size_t i = 0; i < grid->numCells(); ++i ) {
                Cell cell = grid->cell( i );

                if (cell.type() == Cell::Type::TETRAHEDRON){
                    Tetrahedron tetrahedron;
                    tetrahedron.vertices = { cell.index( 0 ), cell.index( 1 ), cell.index( 2 ), cell.index( 3 ) };
                    tetrahedron.splitLevel = 0;

                    gridCells.push_back( tetrahedron );
                } else if (cell.type() == Cell::Type::HEXAHEDRON) {
                    // Decompose cube into five tetrahedra cells
                    std::vector< size_t > vertIndices = { 3, 4, 5, 7, 2, 1, 3, 5, 5, 7, 1, 3, 0, 1, 3, 7, 6, 1, 7, 5 };
                    for ( size_t j = 0; j < vertIndices.size() / 4; ++j ) {
                        Tetrahedron tetrahedron;
                        tetrahedron.splitLevel = 0;
                        tetrahedron.vertices = { 
                            cell.index( vertIndices[ 4 * j + 0 ] ),
                            cell.index( vertIndices[ 4 * j + 1 ] ),
                            cell.index( vertIndices[ 4 * j + 2 ] ),
                            cell.index( vertIndices[ 4 * j + 3 ] ) };

                        gridCells.push_back( tetrahedron );
                    }
                } else {
                    throw std::logic_error( "Grid can only consist of tetrahedron or hexahedron cells!");
                }
            }
            
            debugLog() << "Grid consists of " << gridCells.size() << " tetrahedra cells" << std::endl;
            
            std::vector< PointF< 3 > > isoSurfaceVertices;
            std::vector< uint32_t > isoSurfaceTriangulationIndices;

            // Go through all cells in the grid
            for ( size_t i = 0; i < gridCells.size(); ++i ) {
                auto cellVertices = gridCells[ i ].vertices;
                // Tetrahedron has six edges that need to be added to the hash maps
                for ( uint8_t j = 0; j < 4; ++j ) {
                    for ( uint8_t k = j + 1; k < 4; ++k ) {
                        EdgeKey edgeKey = getEdgeKey( cellVertices[ j ], cellVertices[ k ] );

                        if ( incidentCellsOfEdge.find( edgeKey ) == incidentCellsOfEdge.end() ) {
                            // Edge has no entry yet, initialize a new vector for this entry
                            auto cellList = std::vector< size_t >();
                            incidentCellsOfEdge[ edgeKey ] = cellList;
                        }
                        // Add cell index to vector
                        incidentCellsOfEdge[ edgeKey ].push_back( i );
                    }
                } 
            }


            //==========================================================================================================
            // STEP 2: Form list of diamonds by combining tetrahedron cells that are incident to an edge
            //==========================================================================================================

            for ( size_t i = 0; i < gridCells.size(); ++i ) {
                auto cellVertices = gridCells[ i ].vertices;
                size_t splitLevel = gridCells[ i ].splitLevel;

                // Skip cells that are marked as invalid (deleted cells)
                if ( splitLevel == 255 ) continue;

                for ( uint8_t j = 0; j < 4; ++j ) {
                    for ( uint8_t k = j + 1; k < 4; ++k ) {
                        EdgeKey edgeKey = getEdgeKey( cellVertices[ j ], cellVertices[ k ] );

                        // Skip edges that have already been processed
                        if ( edgeIntersections.find( edgeKey ) != edgeIntersections.end() ) {
                            continue;
                        }

                        // Use Marching Tetrahedra instead of the Marching Diamond procedure if one of three cases occur:
                        // 1) Forming diamonds is disabled by the user options
                        // 2) The cell has reached the split limit
                        // 3) The current edge does not form a diamond with it's incident cells
                        std::vector< size_t > diamondPoints;
                        if ( noDiamonds || splitLevel == maxSplits || !getDiamondPoints( edgeKey, diamondPoints ) ) {
                            // Calculate intersection point based on linear interpolation
                            double s1 = getFieldValue( gridPoints[ edgeKey.first ] );
                            double s2 = getFieldValue( gridPoints[ edgeKey.second ] );
                            if ( ( s1 >= iso && s2 <= iso ) || ( s1 <= iso && s2 >= iso ) ) {
                                // Calculate coefficients to determine the iso point on the line
                                double alpha = std::abs( s1 - iso ) / std::abs( s1 - s2 );
                                auto p1 = gridPoints[ edgeKey.first ];
                                auto p2 = gridPoints[ edgeKey.second ];

                                edgeIntersections[ edgeKey ] = p1 * ( 1 - alpha ) + p2 * alpha;
                            }
                            continue;
                        }

                        // Get corresponding values of diamond vertices from the scalar field
                        uint8_t K = diamondPoints.size() - 2;
                        std::vector< double > scalars( K + 2 );
                        double ringSum = 0;
                        for ( uint8_t j = 0; j < K + 2; ++j ) {
                            scalars[ j ] = getFieldValue( gridPoints[ diamondPoints[ j ] ] );
                            if ( j < K ) {
                                ringSum += scalars[ j ];
                            }
                        }

                        // Calculate constants that depend on K
                        const double C = std::sin( 4 * M_PI / K ) - 2 * std::sin( 2 * M_PI / K );
                        const double D = ( K - 2 ) * std::abs(
                            3 * std::sin( 2 * M_PI / K ) - 4 * std::sin( 4 * M_PI / K ) - std::sin( 6 * M_PI / K ) );

                        // Calculate other constants that depend on the scalar values and the iso value
                        const double L = D * ( scalars[ K ] - iso );
                        const double M = 4 * std::abs( C ) * ( ringSum - K * iso );
                        const double N = D * ( scalars[ K + 1 ] - iso );

                        // Solve roots of cubic polynomial to determine intersection points (0,0,z) in the reference diamond
                        // Cubic polynomial: f(z) = L(2-z)^3 + Mz(2-z)^3 + Nz^3 = 8L + (4M-12L)z + (6L-4M)z^2 + (M+N-L)z^3
                        double coefficients[4] = { 8 * L, 4 * M - 12 * L, 6 * L - 4 * M, M + N - L };
                        double roots[3];
                        int numberOfRoots = solveCubicPolynomial( coefficients, roots );

                        // Only keep roots that inside the intervall [0,2]
                        std::vector< double > rootsInsideDiamond;
                        for ( uint8_t j = 0; j < numberOfRoots; ++j ) {
                            if ( roots[ j ] >= 0 && roots[ j ] <= 2 ) {
                                rootsInsideDiamond.push_back( roots[ j ] );
                            }
                        }

                        if ( rootsInsideDiamond.size() == 0 ) {
                            continue;
                        }

                        double z;
                        if ( rootsInsideDiamond.size() == 1 ) {
                            z = rootsInsideDiamond[ 0 ];
                        } else if ( rootsInsideDiamond.size() == 2 ) {
                            // get midpoint between roots
                            z = ( rootsInsideDiamond[ 0 ] + rootsInsideDiamond[ 1 ] ) / 2;
                        }

                        std::vector< double > barycentricCoords;
                        barycentricCoords.reserve( K + 2 );

                        const double E = 4 * K * std::abs( C ) * z * std::pow( 2 - z, 2 )
                            + D * ( std::pow( 2 - z, 3 ) + std::pow( z, 3 ) );

                        for ( uint8_t k = 0; k < K; ++k ) {
                            barycentricCoords[ k ] = 4 * std::abs( C ) * z * std::pow( 2 - z, 2 ) / E;
                        }
                        barycentricCoords[ K ] = D * std::pow( 2 - z, 3 ) / E;
                        barycentricCoords[ K + 1 ] = D * std::pow( z, 3 ) / E;

                        // Transform barycentric coordinates to point in the field
                        Vector3 zPoint(0, 0, 0);
                        for ( uint8_t j = 0; j < K + 2; ++j ) {
                            zPoint += barycentricCoords[ j ] * gridPoints[ diamondPoints[ j ] ];
                        }

                        if ( rootsInsideDiamond.size() == 1 ) {
                            edgeIntersections[ edgeKey ] = zPoint;
                        } else if ( rootsInsideDiamond.size() == 2 ) {
                            // Found two intersection along the incident line, so the diamon needs to be split
                            // Add split point as a new point to the grid
                            size_t splitPointIndex = gridPoints.size();
                            gridPoints.push_back( zPoint );

                            // Delete cells that are connected to the current incident edge
                            // In order to keep index references to the cells in `gridCells` consistent, the elements are not removed
                            // from the array. Instead the splitLevel attribute is set to 255, which indicates that the cell is invalid
                            for ( size_t oldCellIndices : incidentCellsOfEdge[ edgeKey ] ) {
                                gridCells[ oldCellIndices ].splitLevel = 255;
                            }

                            // Insert new cells
                            size_t cellOffset = gridCells.size();
                            for ( uint8_t n = 0; n < K; ++n ) {
                                Tetrahedron t1, t2;
                                t1.vertices = { splitPointIndex, diamondPoints[ K + 1 ], diamondPoints[ n ], diamondPoints[ ( n + 1 ) % K ] };
                                t2.vertices = { splitPointIndex, diamondPoints[ K ], diamondPoints[ ( n + 1 ) % K ], diamondPoints[ n ] };
                                t1.splitLevel = splitLevel + 1;
                                t2.splitLevel = splitLevel + 1;

                                gridCells.push_back( t1 );
                                gridCells.push_back( t2 );
                            }

                            // Update `incidentCellsOfEdge` map for the newly created edges and cells
                            for ( uint8_t m = 0; m < K + 2; ++m ) {
                                EdgeKey edgeWithSplitPoint = getEdgeKey( splitPointIndex, diamondPoints[ m ] );
                                // Initialize cell list of new edge with empty array
                                incidentCellsOfEdge[ edgeWithSplitPoint ] = std::vector< size_t >();
                                // Add all cells to the list that are incident to the current edge
                                for ( uint8_t n = 0; n < 2 * K; ++n ) {
                                    auto cellVertices = gridCells[ cellOffset + n ].vertices;
                                    // The cell must contain both endpoints of the current edge in the vertices array in order to be incident
                                    if ( std::find( std::begin( cellVertices ), std::end( cellVertices ), splitPointIndex ) == std::end( cellVertices ) ) continue;
                                    if ( std::find( std::begin( cellVertices ), std::end( cellVertices ), diamondPoints[ m ]) == std::end( cellVertices ) ) continue;

                                    incidentCellsOfEdge[ edgeWithSplitPoint ].push_back( cellOffset + n );
                                } 
                            }
                        }
                    }
                }
            }

            //==========================================================================================================
            // STEP 3: Add faces of iso surface to each tetrahedron cell
            //==========================================================================================================

            for ( size_t i = 0; i < gridCells.size(); ++i ) {
                auto cellVertices = gridCells[ i ].vertices;

                // Skip cells that are marked as invalid
                if ( gridCells[ i ].splitLevel == 255 ) continue;
                
                uint8_t numPoints = 0;
                uint32_t offset = (uint32_t)isoSurfaceVertices.size();
                // Tetrahedron has six edges that need to be added to the hash maps
                for ( uint8_t j = 0; j < 4; ++j ) {
                    for ( uint8_t k = j + 1; k < 4; ++k ) {
                        EdgeKey edgeKey = getEdgeKey( cellVertices[ j ], cellVertices[ k ] );

                        if ( edgeIntersections.find( edgeKey ) != edgeIntersections.end() ) {
                            isoSurfaceVertices.push_back( (Vector3F)edgeIntersections[ edgeKey ] );
                            numPoints++;
                        }
                    }
                }

                // Connect points to triangular face
                if ( numPoints == 3 ) {
                    isoSurfaceTriangulationIndices.insert( isoSurfaceTriangulationIndices.end(), {
                        offset, offset + 1, offset + 2 });
                } else if ( numPoints == 4 ) {
                    isoSurfaceTriangulationIndices.insert( isoSurfaceTriangulationIndices.end(), {
                        offset + 0, offset + 1, offset + 2,
                        offset + 2, offset + 1, offset + 3 });
                }
            }

            debugLog() << "Vertices: " << isoSurfaceVertices.size() << " Triangles: " << isoSurfaceTriangulationIndices.size() / 3 << std::endl;

            // Compute bounding sphere once based on the wireframe vertices, this should be the same for the surface vertices
            auto bs = graphics::computeBoundingSphere( isoSurfaceVertices );

            // The GraphicsSystem is needed to create Drawables, which represent the to be rendererd objects.
            auto const& system = graphics::GraphicsSystem::instance();
            std::string resourcePath = PluginRegistrationService::getInstance().getResourcePath( "utils/Graphics" );

            auto surfaceNormals = graphics::computeNormals( isoSurfaceVertices, isoSurfaceTriangulationIndices );
            
            auto isoSurfaceDrawable = system.makePrimitive(
                graphics::PrimitiveConfig{ graphics::RenderPrimitives::TRIANGLES }
                    .vertexBuffer( "position", system.makeBuffer( isoSurfaceVertices ) )
                    .vertexBuffer( "normal", system.makeBuffer( surfaceNormals ) )
                    .indexBuffer( system.makeIndexBuffer( isoSurfaceTriangulationIndices ) )
                    .uniform( "color", Color(1.0f, 0.8f, 0.1f) )
                    .renderOption( graphics::RenderOption::Blend, true )
                    .boundingSphere( bs ),
                system.makeProgramFromFiles( resourcePath + "shader/surface/phong/singleColor/vertex.glsl",
                                             resourcePath + "shader/surface/phong/singleColor/fragment.glsl" ) );

            setGraphics( "Iso surface", isoSurfaceDrawable );
        }
    };
 
    AlgorithmRegister< MarchingDiamondsAlgorithm > dummy( "Hauptaufgabe/1.2.1 Marching Diamonds", "" );
} // namespace

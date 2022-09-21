using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using Unity.Mathematics;

// Adapted from the C++ CEAlign code v0.9: https://pymolwiki.org/index.php/Cealign_plugin
// 
// BSD 3-Clause License

// Copyright (c) 2022, Nanome
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

namespace Nanome
{
    public class CEAlign
    {
        private const int MAX_KEPT = 20;

        private struct Path
        {
            public int2[] P;
            public Path(int len)
            {
                P = new int2[len];
            }
        }

        private struct PathCache
        {
            public Path[] Paths;
            public PathCache(int len)
            {
                Paths = new Path[len];
            }
        }


        public float4x4 Run(float3[] prot1, int lenA, float3[] prot2,
                        int lenB, out int resLen, out float resRMSD, float d0,
                        float d1, int windowSize, int gapMax)
        {

            int smaller;
            int bufferSize;
            PathCache paths = default;

            smaller = lenA < lenB ? lenA : lenB;

            // calculate the distance matrix for each protein
            float[] dmA = calcDM(prot1, lenA);
            float[] dmB = calcDM(prot2, lenB);

            // calculate the CE Similarity matrix
            float[] S = calcS(dmA, dmB, lenA, lenB, windowSize);

            // find the best path through the CE Sim. matrix
            paths = (PathCache)findPath(S, dmA, dmB, lenA, lenB, d0, d1, windowSize, gapMax, out bufferSize);

            float bestRMSD;
            int bestLen = 0;
            // Get the optimal superposition here...
            float4x4 result = findBest(prot1, prot2, paths, bufferSize, smaller, windowSize, out bestLen, out bestRMSD);

            resRMSD = (float)bestRMSD;
            resLen = bestLen;

            return result;
        }


        float[] calcDM(float3[] coords, int len)
        {
            float[] result = new float[len * len];

            for (int row = 0; row < len; row++)
            {
                for (int col = 0; col < len; col++)
                {
                    int id = row * len + col;
                    result[id] = math.distance(coords[row], coords[col]);
                }
            }
            return result;
        }

        float[] calcS(float[] d1, float[] d2, int lenA, int lenB, int wSize)
        {
            int i;
            float winSize = (float)wSize;
            // initialize the 2D similarity matrix
            float[] S = new float[lenA * lenB];

            float sumSize = (winSize - 1.0f) * (winSize - 2.0f) / 2.0f;
            //
            // This is where the magic of CE comes out.  In the similarity matrix,
            // for each i and j, the value of ceSIM[i][j] is how well the residues
            // i - i+winSize in protein A, match to residues j - j+winSize in protein
            // B.  A value of 0 means absolute match; a value >> 1 means bad match.
            //
            int iA, iB, row, col;
            for (iA = 0; iA < lenA; iA++)
            {
                for (iB = 0; iB < lenB; iB++)
                {
                    int id = iA * lenB + iB;
                    S[id] = -1.0f;

                    if (iA > lenA - wSize || iB > lenB - wSize)
                        continue;

                    float score = 0.0f;

                    //
                    // We always skip the calculation of the distance from THIS
                    // residue, to the next residue.  This is a time-saving heur-
                    // istic decision.  Almost all alpha carbon bonds of neighboring
                    // residues is 3.8 Angstroms.  Due to entropy, S = -k ln pi * pi,
                    // this tell us nothing, so it doesn't help so ignore it.
                    //
                    for (row = 0; row < wSize - 2; row++)
                    {
                        for (col = row + 2; col < wSize; col++)
                        {
                            int id1 = (iA + row) * lenA + (iA + col);
                            int id2 = (iB + row) * lenB + (iB + col);
                            score += math.abs(d1[id1] - d2[id2]);
                        }
                    }

                    S[id] = score / sumSize;
                }
            }
            return S;
        }


        PathCache findPath(float[] S, float[] dA, float[] dB, int lenA, int lenB, float D0, float D1, int winSize, int gapMax, out int bufferSize)
        {
            // the best Path's score
            float bestPathScore = 1e6f;
            int bestPathLength = 0;

            // length of longest possible alignment
            int smallerSize = (lenA < lenB) ? lenA : lenB;
            int winSum = (winSize - 1) * (winSize - 2) / 2;

            Path bestPath = new Path(smallerSize);

            // index variable for below
            for (int i = 0; i < smallerSize; i++)
            {
                bestPath.P[i].x = -1;
                bestPath.P[i].y = -1;
            }

            //======================================================================
            // for storing the best 20 paths
            bufferSize = 0;
            int bufferIndex = 0;
            int[] lenBuffer = new int[MAX_KEPT];
            float[] scoreBuffer = new float[MAX_KEPT];
            PathCache pathBuffer = new PathCache(MAX_KEPT);

            for (int i = 0; i < MAX_KEPT; i++)
            {
                // initialize the paths
                scoreBuffer[i] = 1e6f;
                lenBuffer[i] = 0;
                pathBuffer.Paths[i] = new Path(0);
            }

            // winCache
            // this array stores a list of residues seen.  We use it to calculate the
            // total score of a path from 1..M and then add it to M+1..N.
            int[] winCache = new int[smallerSize];
            for (int i = 0; i < smallerSize; i++)
            {
                winCache[i] = (i + 1) * i * winSize / 2 + (i + 1) * winSum;
            }

            // allScoreBuffer
            // this 2D array keeps track of all partial gapped scores
            int scoreBufferM = gapMax * 2 + 1;
            float[] allScoreBuffer = new float[smallerSize * scoreBufferM];

            // initialize the ASB
            for (int i = 0; i < allScoreBuffer.Length; i++)
            {
                allScoreBuffer[i] = 1e6f;
            }

            int[] tIndex = new int[smallerSize];
            int gapBestIndex = -1;

            //======================================================================
            // Start the search through the CE matrix.
            //
            for (int iA = 0; iA < lenA; iA++)
            {
                if (iA > lenA - winSize * (bestPathLength - 1))
                    break;

                for (int iB = 0; iB < lenB; iB++)
                {
                    int id = iA * lenB + iB;
                    if (S[id] >= D0)
                        continue;

                    if (S[id] == -1.0f)
                        continue;

                    if (iB > lenB - winSize * (bestPathLength - 1))
                        break;

                    //
                    // Restart curPath here.
                    //
                    Path curPath = new Path(smallerSize);

                    for (int a = 0; a < smallerSize; a++)
                    {
                        curPath.P[a].x = -1;
                        curPath.P[a].y = -1;
                    }
                    curPath.P[0].x = iA;
                    curPath.P[0].y = iB;
                    int curPathLength = 1;
                    tIndex[0] = 0;
                    float curTotalScore = 0.0f;

                    //
                    // Check all possible paths starting from iA, iB
                    //
                    bool done = false;
                    while (!done)
                    {
                        float gapBestScore = 1e6f;
                        gapBestIndex = -1;

                        // Check all possible gaps [1..gapMax] from here
                        for (int g = 0; g < (gapMax * 2) + 1; g++)
                        {
                            int jA = curPath.P[curPathLength - 1].x + winSize;
                            int jB = curPath.P[curPathLength - 1].y + winSize;

                            if ((g + 1) % 2 == 0)
                            {
                                jA += (g + 1) / 2;
                            }
                            else
                            { // ( g odd )
                                jB += (g + 1) / 2;
                            }

                            //
                            // Following are three heuristics to ensure high quality
                            // long paths and make sure we don't run over the end of
                            // the S, matrix.

                            // 1st: If jA and jB are at the end of the matrix
                            if (jA > lenA - winSize || jB > lenB - winSize)
                            {
                                // FIXME, was: jA > lenA-winSize-1 || jB > lenB-winSize-1
                                continue;
                            }
                            // 2nd: If this gapped octapeptide is bad, ignore it.
                            if (S[jA * lenB + jB] > D0)
                                continue;
                            // 3rd: if too close to end, ignore it.
                            if (S[jA * lenB + jB] == -1.0f)
                                continue;

                            float curScore = 0.0f;
                            int s;
                            for (s = 0; s < curPathLength; s++)
                            {
                                int id1 = curPath.P[s].x * lenA + jA;
                                int id2 = curPath.P[s].y * lenB + jB;
                                curScore += math.abs(dA[id1] - dB[id2]);

                                int id3 = (curPath.P[s].x + (winSize - 1)) * lenA + (jA + (winSize - 1));
                                int id4 = (curPath.P[s].y + (winSize - 1)) * lenB + (jB + (winSize - 1));
                                curScore += math.abs(dA[id3] - dB[id4]);

                                for (int k = 1; k < winSize - 1; k++)
                                {
                                    int id5 = (curPath.P[s].x + k) * lenA + (jA + (winSize - 1) - k);
                                    int id6 = (curPath.P[s].y + k) * lenB + (jB + (winSize - 1) - k);
                                    curScore += math.abs(dA[id5] - dB[id6]);
                                }
                            }

                            curScore /= (float)(winSize * curPathLength);

                            if (curScore >= D1)
                            {
                                continue;
                            }

                            // store GAPPED best
                            if (curScore < gapBestScore)
                            {
                                curPath.P[curPathLength].x = jA;
                                curPath.P[curPathLength].y = jB;
                                gapBestScore = curScore;
                                gapBestIndex = g;
                                allScoreBuffer[(curPathLength - 1) * scoreBufferM + g] = curScore;
                            }
                        } /// ROF -- END GAP SEARCHING

                        //
                        // DONE GAPPING:
                        //

                        // calculate curTotalScore
                        curTotalScore = 0.0f;
                        int jGap, gA, gB;
                        float score1 = 0.0f, score2 = 0.0f;

                        if (gapBestIndex != -1)
                        {
                            jGap = (gapBestIndex + 1) / 2;
                            if ((gapBestIndex + 1) % 2 == 0)
                            {
                                gA = curPath.P[curPathLength - 1].x + winSize + jGap;
                                gB = curPath.P[curPathLength - 1].y + winSize;
                            }
                            else
                            {
                                gA = curPath.P[curPathLength - 1].x + winSize;
                                gB = curPath.P[curPathLength - 1].y + winSize + jGap;
                            }

                            int idgAgB = gA * lenB + gB;
                            int idScore1 = (curPathLength - 1) * scoreBufferM + gapBestIndex;
                            // perfect
                            score1 = (allScoreBuffer[idScore1] * winSize * curPathLength + S[idgAgB] * winSum) / (winSize * curPathLength + winSum);

                            int idScore2 = (curPathLength - 2) * scoreBufferM + tIndex[curPathLength - 1];
                            // perfect
                            score2 = ((curPathLength > 1 ? (allScoreBuffer[idScore2])
                                                        : S[id]) *
                                        winCache[curPathLength - 1] +
                                    score1 * (winCache[curPathLength] - winCache[curPathLength - 1])) /
                                    winCache[curPathLength];

                            curTotalScore = score2;
                            // heuristic -- path is getting sloppy, stop looking
                            if (curTotalScore > D1)
                            {
                                done = true;
                                gapBestIndex = -1;
                                break;
                            }
                            else
                            {
                                allScoreBuffer[(curPathLength - 1) * scoreBufferM + gapBestIndex] = curTotalScore;
                                tIndex[curPathLength] = gapBestIndex;
                                curPathLength++;
                            }
                        }
                        else
                        {
                            // if here, then there was no good gapped path
                            // so quit and restart from iA, iB+1
                            done = true;
                            curPathLength--;
                            break;
                        }

                        //
                        // test this gapped path against the best seen
                        // starting from iA, iB
                        //

                        // if our currently best gapped path from iA and iB is LONGER
                        // than the current best; or, it's equal length and the score's
                        // better, keep the new path.
                        if (curPathLength > bestPathLength ||
                            (curPathLength == bestPathLength && curTotalScore < bestPathScore))
                        {
                            bestPathLength = curPathLength;
                            bestPathScore = curTotalScore;
                            // deep copy curPath
                            Path tempPath = new Path(smallerSize);

                            for (int a = 0; a < smallerSize; a++)
                            {
                                tempPath.P[a].x = curPath.P[a].x;
                                tempPath.P[a].y = curPath.P[a].y;
                            }

                            bestPath = tempPath;
                        }
                    } /// END WHILE

                    //
                    // At this point, we've found the best path starting at iA, iB.
                    //
                    if (bestPathLength > lenBuffer[bufferIndex] ||
                        (bestPathLength == lenBuffer[bufferIndex] && bestPathScore < scoreBuffer[bufferIndex]))
                    {

                        // we're going to add an entry to the ring-buffer.
                        // Adjust maxSize values and curIndex accordingly.
                        bufferIndex = (bufferIndex == MAX_KEPT - 1) ? 0 : bufferIndex + 1;
                        bufferSize = (bufferSize < MAX_KEPT) ? (bufferSize) + 1 : MAX_KEPT;
                        Path pathCopy = new Path(smallerSize);

                        for (int a = 0; a < smallerSize; a++)
                        {
                            pathCopy.P[a].x = bestPath.P[a].x;
                            pathCopy.P[a].y = bestPath.P[a].y;
                        }

                        if (bufferIndex == 0 && bufferSize == MAX_KEPT)
                        {
                            pathBuffer.Paths[MAX_KEPT - 1] = pathCopy;
                            scoreBuffer[MAX_KEPT - 1] = bestPathScore;
                            lenBuffer[MAX_KEPT - 1] = bestPathLength;
                        }
                        else
                        {
                            pathBuffer.Paths[bufferIndex - 1] = pathCopy;
                            scoreBuffer[bufferIndex - 1] = bestPathScore;
                            lenBuffer[bufferIndex - 1] = bestPathLength;
                        }
                    }
                } // ROF -- end for iB
            }     // ROF -- end for iA

            return pathBuffer;
        }


        float4x4 findBest(float3[] coordsA, float3[] coordsB, PathCache paths, int bufferSize, int smallerSize, int winSize, out int bestLen, out float bestRMSD)
        {
            // keep the best values
            bestRMSD = 1e6f;
            bestLen = 0;

            Matrix<float> bestU = Matrix<float>.Build.Dense(3, 3);
            float3 bestCOM1 = float3.zero;
            float3 bestCOM2 = float3.zero;
            // TA2<float> bestU;
            // TA1<float> bestCOM1, bestCOM2;
            int bestO = -1;

            // loop through the buffer
            for (int o = 0; o < bufferSize; o++)
            {

                // grab the current path
                float3[] c1 = new float3[smallerSize];
                float3[] c2 = new float3[smallerSize];

                for (int a = 0; a < smallerSize; a++)
                {
                    c1[a] = float3.zero;
                    c2[a] = float3.zero;
                }
                int curLen = 0;

                int it = 0;
                for (int j = 0; j < smallerSize; j++)
                {

                    Path patho = paths.Paths[o];
                    // rebuild the coordinate lists for this path
                    if (patho.P[j].x != -1)
                    {
                        for (int k = 0; k < winSize; k++)
                        {

                            float3 t1 = coordsA[patho.P[j].x + k].xyz;
                            float3 t2 = coordsB[patho.P[j].y + k].xyz;

                            c1[it] = t1;
                            c2[it] = t2;

                            it++;
                        }
                    }
                    else
                    {
                        curLen = it;
                        break;
                    }
                }

                //
                // For convenience, let there be M points of N dimensions
                //
                int m = curLen;

                //==========================================================================
                //
                // Superpose the two proteins
                //
                //==========================================================================

                // centers of mass for c1 and c2
                float3 c1COM = float3.zero;
                float3 c2COM = float3.zero;

                // Calc CsOM
                for (int i = 0; i < m; i++)
                {
                    c1COM += c1[i];
                    c2COM += c2[i];
                }
                c1COM /= m;
                c2COM /= m;


                // Move the two vectors to the origin
                for (int i = 0; i < m; i++)
                {
                    c1[i] -= c1COM;
                    c2[i] -= c2COM;
                }

                //==========================================================================
                //
                // Calculate U and RMSD.  This is broken down to the super-silly-easy
                // math of: U = Wt * V, where Wt and V are NxN matrices from the SVD of
                // R, the correlation matrix between the two origin-based vector sets.
                //
                //==========================================================================

                // Calculate the initial residual, E0
                // E0 = sum( Yn*Yn + Xn*Xn ) -- sum of squares
                float E0 = 0.0f;
                for (int i = 0; i < m; i++)
                {
                    E0 += (c1[i].x * c1[i].x) + (c2[i].x * c2[i].x);
                    E0 += (c1[i].y * c1[i].y) + (c2[i].y * c2[i].y);
                    E0 += (c1[i].z * c1[i].z) + (c2[i].z * c2[i].z);
                }

                //
                // SVD is the SVD of the correlation matrix Xt*Y
                // R = c2' * c1 = W * S * Vt

                var c1Matrix = Matrix<float>.Build.Dense(c1.Length, 3);
                var c2Matrix = Matrix<float>.Build.Dense(c2.Length, 3);

                for (int i = 0; i < c1.Length; i++)
                {
                    c1Matrix[i, 0] = c1[i].x;
                    c1Matrix[i, 1] = c1[i].y;
                    c1Matrix[i, 2] = c1[i].z;
                }

                for (int i = 0; i < c2.Length; i++)
                {
                    c2Matrix[i, 0] = c2[i].x;
                    c2Matrix[i, 1] = c2[i].y;
                    c2Matrix[i, 2] = c2[i].z;
                }

                var toSVD = c2Matrix.TransposeThisAndMultiply(c1Matrix);

                var svdResult = toSVD.Svd(true);
                var VT = svdResult.VT;
                var U = svdResult.U;
                Vector<float> sigmas = svdResult.S;

                //
                // Check any reflections before rotation of the points;
                // if det(W)*det(V) == -1 then we just reflect
                // the principal axis corresponding to the smallest eigenvalue by -1
                //
                LU<float> LU_VT = VT.LU();
                LU<float> LU_U = U.LU();

                if (LU_U.Determinant * LU_VT.Determinant < 0.0f)
                {
                    // std::cout << "_________REFLECTION_________" << std::endl;

                    // revese the smallest axes and last sigma

                    for (int i = 0; i < 3; i++)
                    {
                        U[2, i] = -U[2, i];
                    }

                    sigmas[2] = -sigmas[2];
                }

                // calculate the rotation matrix, U.
                // U = W * Vt
                var UVT = U.Multiply(VT);

                // Now calculate the RMSD
                float3 sigmas3 = new float3(sigmas[0], sigmas[1], sigmas[2]);
                float sig = math.csum(sigmas3);

                float curRMSD = math.sqrt(math.abs((E0 - 2 * sig) / (float)m));

                //
                // Save the best
                //
                if (curRMSD < bestRMSD || (curRMSD == bestRMSD && c1.Length > bestLen))
                {
                    UVT.CopyTo(bestU);
                    bestRMSD = curRMSD;
                    bestCOM1 = new float3(c1COM);
                    bestCOM2 = new float3(c2COM);
                    bestLen = curLen;
                    bestO = o;
                }
            }

            if (bestRMSD == 1e6f)
            {
                // Could not find a good rotation matrix
                return float4x4.identity;
            }

            float4x4 rVal = new float4x4(
                bestU[0, 0],
                bestU[1, 0],
                bestU[2, 0],
                bestCOM1.x,
                bestU[0, 1],
                bestU[1, 1],
                bestU[2, 1],
                bestCOM1.y,
                bestU[0, 2],
                bestU[1, 2],
                bestU[2, 2],
                bestCOM1.z,
                -bestCOM2[0],
                -bestCOM2[1],
                -bestCOM2[2],
                1.0f
            );

            return rVal;
        }
    }
}

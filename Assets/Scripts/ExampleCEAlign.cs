using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Unity.Mathematics;
using UnityEngine;

namespace Nanome
{
    public class ExampleCEAlign : MonoBehaviour
    {

        void Start()
        {
            string path = Path.Combine(Application.streamingAssetsPath, "Coords_3EI0_3EAM.txt");
            (float3[] coords1, float3[] coords2) = readCoordinates(path);

            CEAlignResult result = CEAlign.Run(coords1, coords2);

            Debug.Log("RMSD = " + result.RMSD);

        }

        private (float3[], float3[]) readCoordinates(string path)
        {
            List<float3> coords1 = new List<float3>();
            List<float3> coords2 = new List<float3>();

            string content = File.ReadAllText(path);
            string[] lines = content.Split(new string[] { "\n", "\r" }, System.StringSplitOptions.RemoveEmptyEntries);

            bool secondSet = false;

            foreach (string l in lines)
            {
                if (l.Length < 3) continue;

                if (l.Contains("==="))
                {
                    secondSet = true;
                    continue;
                }
                string[] tokens = l.Split(new string[] { " ", "\t" }, System.StringSplitOptions.RemoveEmptyEntries);

                float x = float.Parse(tokens[0], CultureInfo.InvariantCulture);
                float y = float.Parse(tokens[1], CultureInfo.InvariantCulture);
                float z = float.Parse(tokens[2], CultureInfo.InvariantCulture);
                float3 position = new float3(x, y, z);
                if (!secondSet)
                {
                    coords1.Add(position);
                }
                else
                {
                    coords2.Add(position);
                }
            }


            return (coords1.ToArray(), coords2.ToArray());
        }
    }
}
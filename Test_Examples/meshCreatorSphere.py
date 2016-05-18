from steps.utilities.meshio import importAbaqus, stetmesh, saveMesh

MyMesh = importAbaqus('Sphere', 1e-7)[0]
print type(MyMesh)
saveMesh('testMeshSphere', MyMesh)

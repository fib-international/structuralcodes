import _utils2


y = (0, 0, 0.35 , 0.35)
z = (0, -0.45, -0.45, 0)

area = _utils2.marin_integration('yz', y, z, 0, 0)
first_moment_yy = _utils2.marin_integration('yz', y, z, 1, 0)
first_moment_zz = _utils2.marin_integration('yz', y, z, 0, 1)

print(f'Area of a rectangular section 0.35x0.45 = {area:.4f}')
print(f'First area of moment yy of a rectangular section 0.35x0.45 = {first_moment_yy:.3f}')
print(f'Centroid y = {first_moment_yy/area:.3f}')
print(f'First area of moment zz of a rectangular section 0.35x0.45 = {first_moment_zz:.3f}')
print(f'Centroid y = {first_moment_zz/area:.3f}')
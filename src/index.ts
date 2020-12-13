/*
 * Returns a list of 2-D vertex coordinates that will create
 * a rectangle, centered at the specified position.
 */
export function generateRectangleVertices(x, y, w, h) {
  return [
    x - w / 2,
    y - h / 2,
    x + w / 2,
    y - h / 2,
    x + w / 2,
    y + h / 2,

    x - w / 2,
    y - h / 2,
    x + w / 2,
    y + h / 2,
    x - w / 2,
    y + h / 2,
  ];
}

export function getVerts(width, length, height) {
  return [
    // Front
    [-width, length, height], // v0
    [width, length, height], // v1
    [width, length, -height], // v2
    [-width, length, -height], // v3

    // Back
    [width, -length, height], // v4
    [-width, -length, height], // v5
    [-width, -length, -height], // v6
    [width, -length, -height], // v7

    // Left
    [width, length, height], // v1
    [width, -length, height], // v4
    [width, -length, -height], // v7
    [width, length, -height], // v2

    // Right
    [-width, -length, height], // v5
    [-width, length, height], // v0
    [-width, length, -height], // v3
    [-width, -length, -height], // v6

    // Top
    [width, length, height], // v1
    [-width, length, height], // v0
    [-width, -length, height], // v5
    [width, -length, height], // v4

    // Bottom
    [width, -length, -height], // v7
    [-width, -length, -height], // v6
    [-width, length, -height], // v3
    [width, length, -height], // v2
  ];
}

export function generateRectangleTexcoords() {
  return [0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1];
}

export function matrixIdentity3x3() {
  const arr = new Float32Array(9);
  arr[0] = 1;
  arr[4] = 1;
  arr[8] = 1;
  return arr;
}

export function matrixCopy3x3(src) {
  return [
    src[0],
    src[1],
    src[2],
    src[3],
    src[4],
    src[5],
    src[6],
    src[7],
    src[8],
  ];
}

export function matrixMultiply3x3(...args) {
  if (args.length === 0) {
    throw new Error('At least two matrices must be provided.');
  }
  if (args.length === 1) {
    return args[0];
  }
  let rv = matrixCopy3x3(args[0]);
  for (let i = 1; i < args.length; ++i) {
    const a = rv;
    const b = args[i];
    rv = [
      a[0] * b[0] + a[1] * b[3] + a[2] * b[6],
      a[0] * b[1] + a[1] * b[4] + a[2] * b[7],
      a[0] * b[2] + a[1] * b[5] + a[2] * b[8],

      a[3] * b[0] + a[4] * b[3] + a[5] * b[6],
      a[3] * b[1] + a[4] * b[4] + a[5] * b[7],
      a[3] * b[2] + a[4] * b[5] + a[5] * b[8],

      a[6] * b[0] + a[7] * b[3] + a[8] * b[6],
      a[6] * b[1] + a[7] * b[4] + a[8] * b[7],
      a[6] * b[2] + a[7] * b[5] + a[8] * b[8],
    ];
  }

  return rv;
}

export function matrixMultiply3x3I(dest, a, b) {
  return matrixSet3x3(
      dest,
      a[0] * b[0] + a[1] * b[3] + a[2] * b[6],
      a[0] * b[1] + a[1] * b[4] + a[2] * b[7],
      a[0] * b[2] + a[1] * b[5] + a[2] * b[8],
      a[3] * b[0] + a[4] * b[3] + a[5] * b[6],
      a[3] * b[1] + a[4] * b[4] + a[5] * b[7],
      a[3] * b[2] + a[4] * b[5] + a[5] * b[8],
      a[6] * b[0] + a[7] * b[3] + a[8] * b[6],
      a[6] * b[1] + a[7] * b[4] + a[8] * b[7],
      a[6] * b[2] + a[7] * b[5] + a[8] * b[8],
  );
}

export function matrixTransform2D(world, x, y) {
  return [
    world[0] * x + world[1] * y + world[2],
    world[3] * x + world[4] * y + world[5],
  ];
}

export function makeTranslation3x3(tx, ty) {
  return makeTranslation3x3I(matrixIdentity3x3(), tx, ty);
}

export function makeTranslation3x3I(dest, tx, ty) {
  return matrixSet3x3(dest, 1, 0, 0, 0, 1, 0, tx, ty, 1);
}

export function makeRotation3x3(angleInRadians) {
  const c = Math.cos(angleInRadians);
  const s = Math.sin(angleInRadians);
  return [c, -s, 0, s, c, 0, 0, 0, 1];
}

export function makeScale3x3(sx, sy) {
  if (arguments.length === 1) {
    sy = sx;
  }
  return makeScale3x3I(matrixIdentity3x3(), sx, sy);
}

export function matrixSet3x3(dest, a1, a2, a3, a4, a5, a6, a7, a8, a9) {
  dest[0] = a1;
  dest[1] = a2;
  dest[2] = a3;
  dest[3] = a4;
  dest[4] = a5;
  dest[5] = a6;
  dest[6] = a7;
  dest[7] = a8;
  dest[8] = a9;
  return dest;
}

export function makeScale3x3I(dest, sx, sy) {
  if (arguments.length === 2) {
    sy = sx;
  }
  return matrixSet3x3(dest, sx, 0, 0, 0, sy, 0, 0, 0, 1);
}

// http://stackoverflow.com/questions/983999/simple-3x3-matrix-inverse-code-c
export function makeInverse3x3(input) {
  const m = function(col, row) {
    return input[row * 3 + col];
  };
  // computes the inverse of a matrix m
  const det =
    m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
    m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
    m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));

  const invdet = 1 / det;

  return [
    (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invdet,
    (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invdet,
    (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invdet,
    (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invdet,
    (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invdet,
    (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invdet,
    (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invdet,
    (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invdet,
    (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invdet,
  ];
}

export function midPoint(x1, y1, x2, y2) {
  return [x1 + (x2 - x1) * 0.5, y1 + (y2 - y1) * 0.5];
}

let vflip = false;
export function getVFlip() {
  return vflip;
}

export function setVFlip(value) {
  vflip = !!value;
}

export function flipVFlip() {
  vflip = !value;
}

export function make2DProjection(width, height, flipVertical) {
  if (flipVertical === undefined) {
    flipVertical = getVFlip();
  }
  flipVertical = flipVertical === true;
  // console.log("Making 2D projection (flipVertical=" + flipVertical + ")");
  flipVertical = flipVertical ? -1 : 1;
  return [
    2 / width,
    0,
    0,
    0,
    -2 / (flipVertical * height),
    0,
    -1,
    flipVertical,
    1,
  ];
}

export function subtractVectors3D(a, b) {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

export function normalize3D(v) {
  const length = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  // make sure we don't divide by 0.
  if (length > 0.00001) {
    return [v[0] / length, v[1] / length, v[2] / length];
  } else {
    return [0, 0, 0];
  }
}

export function cross3D(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

export function makePerspective(fieldOfViewInRadians, aspect, near, far) {
  const f = Math.tan(Math.PI * 0.5 - 0.5 * fieldOfViewInRadians);
  const rangeInv = 1.0 / (near - far);

  return [
    f / aspect,
    0,
    0,
    0,
    0,
    f,
    0,
    0,
    0,
    0,
    (near + far) * rangeInv,
    -1,
    0,
    0,
    near * far * rangeInv * 2,
    0,
  ];
}

export function makeTranslation4x4(tx, ty, tz) {
  return [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, tx, ty, tz, 1];
}

export function makeXRotation(angleInRadians) {
  const c = Math.cos(angleInRadians);
  const s = Math.sin(angleInRadians);

  return [1, 0, 0, 0, 0, c, s, 0, 0, -s, c, 0, 0, 0, 0, 1];
}

export function makeYRotation(angleInRadians) {
  const c = Math.cos(angleInRadians);
  const s = Math.sin(angleInRadians);

  return [c, 0, -s, 0, 0, 1, 0, 0, s, 0, c, 0, 0, 0, 0, 1];
}

export function makeZRotation(angleInRadians) {
  const c = Math.cos(angleInRadians);
  const s = Math.sin(angleInRadians);
  return [c, s, 0, 0, -s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
}

export function makeScale4x4(sx, sy, sz) {
  return [sx, 0, 0, 0, 0, sy, 0, 0, 0, 0, sz, 0, 0, 0, 0, 1];
}

export function matrixMultiply4x4(a, b) {
  const a00 = a[0 * 4 + 0];
  const a01 = a[0 * 4 + 1];
  const a02 = a[0 * 4 + 2];
  const a03 = a[0 * 4 + 3];
  const a10 = a[1 * 4 + 0];
  const a11 = a[1 * 4 + 1];
  const a12 = a[1 * 4 + 2];
  const a13 = a[1 * 4 + 3];
  const a20 = a[2 * 4 + 0];
  const a21 = a[2 * 4 + 1];
  const a22 = a[2 * 4 + 2];
  const a23 = a[2 * 4 + 3];
  const a30 = a[3 * 4 + 0];
  const a31 = a[3 * 4 + 1];
  const a32 = a[3 * 4 + 2];
  const a33 = a[3 * 4 + 3];
  const b00 = b[0 * 4 + 0];
  const b01 = b[0 * 4 + 1];
  const b02 = b[0 * 4 + 2];
  const b03 = b[0 * 4 + 3];
  const b10 = b[1 * 4 + 0];
  const b11 = b[1 * 4 + 1];
  const b12 = b[1 * 4 + 2];
  const b13 = b[1 * 4 + 3];
  const b20 = b[2 * 4 + 0];
  const b21 = b[2 * 4 + 1];
  const b22 = b[2 * 4 + 2];
  const b23 = b[2 * 4 + 3];
  const b30 = b[3 * 4 + 0];
  const b31 = b[3 * 4 + 1];
  const b32 = b[3 * 4 + 2];
  const b33 = b[3 * 4 + 3];
  return [
    a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30,
    a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31,
    a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32,
    a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33,
    a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30,
    a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31,
    a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32,
    a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33,
    a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30,
    a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31,
    a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32,
    a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33,
    a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30,
    a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31,
    a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32,
    a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33,
  ];
}

export function makeInverse4x4(m) {
  const m00 = m[0 * 4 + 0];
  const m01 = m[0 * 4 + 1];
  const m02 = m[0 * 4 + 2];
  const m03 = m[0 * 4 + 3];
  const m10 = m[1 * 4 + 0];
  const m11 = m[1 * 4 + 1];
  const m12 = m[1 * 4 + 2];
  const m13 = m[1 * 4 + 3];
  const m20 = m[2 * 4 + 0];
  const m21 = m[2 * 4 + 1];
  const m22 = m[2 * 4 + 2];
  const m23 = m[2 * 4 + 3];
  const m30 = m[3 * 4 + 0];
  const m31 = m[3 * 4 + 1];
  const m32 = m[3 * 4 + 2];
  const m33 = m[3 * 4 + 3];
  const tmp0 = m22 * m33;
  const tmp1 = m32 * m23;
  const tmp2 = m12 * m33;
  const tmp3 = m32 * m13;
  const tmp4 = m12 * m23;
  const tmp5 = m22 * m13;
  const tmp6 = m02 * m33;
  const tmp7 = m32 * m03;
  const tmp8 = m02 * m23;
  const tmp9 = m22 * m03;
  const tmp10 = m02 * m13;
  const tmp11 = m12 * m03;
  const tmp12 = m20 * m31;
  const tmp13 = m30 * m21;
  const tmp14 = m10 * m31;
  const tmp15 = m30 * m11;
  const tmp16 = m10 * m21;
  const tmp17 = m20 * m11;
  const tmp18 = m00 * m31;
  const tmp19 = m30 * m01;
  const tmp20 = m00 * m21;
  const tmp21 = m20 * m01;
  const tmp22 = m00 * m11;
  const tmp23 = m10 * m01;

  const t0 =
    tmp0 * m11 +
    tmp3 * m21 +
    tmp4 * m31 -
    (tmp1 * m11 + tmp2 * m21 + tmp5 * m31);
  const t1 =
    tmp1 * m01 +
    tmp6 * m21 +
    tmp9 * m31 -
    (tmp0 * m01 + tmp7 * m21 + tmp8 * m31);
  const t2 =
    tmp2 * m01 +
    tmp7 * m11 +
    tmp10 * m31 -
    (tmp3 * m01 + tmp6 * m11 + tmp11 * m31);
  const t3 =
    tmp5 * m01 +
    tmp8 * m11 +
    tmp11 * m21 -
    (tmp4 * m01 + tmp9 * m11 + tmp10 * m21);

  const d = 1.0 / (m00 * t0 + m10 * t1 + m20 * t2 + m30 * t3);

  return [
    d * t0,
    d * t1,
    d * t2,
    d * t3,
    d *
      (tmp1 * m10 +
        tmp2 * m20 +
        tmp5 * m30 -
        (tmp0 * m10 + tmp3 * m20 + tmp4 * m30)),
    d *
      (tmp0 * m00 +
        tmp7 * m20 +
        tmp8 * m30 -
        (tmp1 * m00 + tmp6 * m20 + tmp9 * m30)),
    d *
      (tmp3 * m00 +
        tmp6 * m10 +
        tmp11 * m30 -
        (tmp2 * m00 + tmp7 * m10 + tmp10 * m30)),
    d *
      (tmp4 * m00 +
        tmp9 * m10 +
        tmp10 * m20 -
        (tmp5 * m00 + tmp8 * m10 + tmp11 * m20)),
    d *
      (tmp12 * m13 +
        tmp15 * m23 +
        tmp16 * m33 -
        (tmp13 * m13 + tmp14 * m23 + tmp17 * m33)),
    d *
      (tmp13 * m03 +
        tmp18 * m23 +
        tmp21 * m33 -
        (tmp12 * m03 + tmp19 * m23 + tmp20 * m33)),
    d *
      (tmp14 * m03 +
        tmp19 * m13 +
        tmp22 * m33 -
        (tmp15 * m03 + tmp18 * m13 + tmp23 * m33)),
    d *
      (tmp17 * m03 +
        tmp20 * m13 +
        tmp23 * m23 -
        (tmp16 * m03 + tmp21 * m13 + tmp22 * m23)),
    d *
      (tmp14 * m22 +
        tmp17 * m32 +
        tmp13 * m12 -
        (tmp16 * m32 + tmp12 * m12 + tmp15 * m22)),
    d *
      (tmp20 * m32 +
        tmp12 * m02 +
        tmp19 * m22 -
        (tmp18 * m22 + tmp21 * m32 + tmp13 * m02)),
    d *
      (tmp18 * m12 +
        tmp23 * m32 +
        tmp15 * m02 -
        (tmp22 * m32 + tmp14 * m02 + tmp19 * m12)),
    d *
      (tmp22 * m22 +
        tmp16 * m02 +
        tmp21 * m12 -
        (tmp20 * m12 + tmp23 * m22 + tmp17 * m02)),
  ];
}

export function matrixVectorMultiply4x4(v, m) {
  const dst = [];
  for (let i = 0; i < 4; ++i) {
    dst[i] = 0.0;
    for (let j = 0; j < 4; ++j) dst[i] += v[j] * m[j * 4 + i];
  }
  return dst;
}

/*
 * Returns a 4x4 matrix that, positioned from the camera position,
 * looks at the target, a position in 3-space, angled using the
 * up vector.
 */
export function makeLookAt(cameraPosition, target, up) {
  const zAxis = normalize3D(subtractVectors3D(cameraPosition, target));
  const xAxis = cross3D(up, zAxis);
  const yAxis = cross3D(zAxis, xAxis);

  return [
    xAxis[0],
    xAxis[1],
    xAxis[2],
    0,
    yAxis[0],
    yAxis[1],
    yAxis[2],
    0,
    zAxis[0],
    zAxis[1],
    zAxis[2],
    0,
    cameraPosition[0],
    cameraPosition[1],
    cameraPosition[2],
    1,
  ];
}

// End methods from webglfundamentals.org



var assert = require("assert");
import { generateRectangleVertices } from "../dist/parsegraph-matrix";

describe("generateRectangleVertices", function () {
  it("works", () => {
    assert.ok(generateRectangleVertices(10, 10, 100, 100));
  });
});

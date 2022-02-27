var assert = require("assert");
import { generateRectangleVertices } from "../src/index";

describe("generateRectangleVertices", function () {
  it("works", () => {
    assert.ok(generateRectangleVertices(10, 10, 100, 100));
  });
});

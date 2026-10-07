import { defineConfig } from "vite";

// Bundle each widget as a single self-contained ES module (plus extracted CSS)
// into the Python package, where anywidget loads it via `_esm` / `_css`.
export default defineConfig({
  build: {
    lib: {
      entry: "src/voter_vignette.ts",
      formats: ["es"],
      fileName: () => "voter_vignette.js",
      cssFileName: "voter_vignette",
    },
    outDir: "../src/valency_anndata/viz/static",
    emptyOutDir: true,
    minify: true,
    // Library mode skips whitespace minification for ES output by default;
    // force it, since this bundle ships inside the Python wheel.
    rolldownOptions: { output: { minify: true } },
  },
});

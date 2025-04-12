import { defineConfig } from 'tsup';

export default defineConfig({
  entry: ['src/index.ts'],
  format: ['cjs', 'esm', 'iife'],
  globalName: 'lamejs2',
  dts: true,
  splitting: false,
  sourcemap: true,
  clean: true,
  minify: true,
  esbuildOptions(options) {
    options.resolveExtensions = ['.ts', '.js'];
    options.platform = 'browser';
    options.bundle = true;
    options.mainFields = ['module', 'main'];
    options.target = 'es2018';
  }
}); 
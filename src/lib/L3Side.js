import { SBMAX_s} from './Encoder.js';

class L3Side {
	/**
	 * max scalefactor band, max(SBMAX_l, SBMAX_s*3, (SBMAX_s-3)*3+8)
	 */
	static SFBMAX = (SBMAX_s * 3);
}

export { L3Side };

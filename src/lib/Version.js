class Version {
    constructor(context) {
        this.context = context;
    }

    // Move all methods to the class definition
    getLameShortVersion() {
        return "3.100";
    }

    getLameVersion() {
        return "3.100";
    }

    getLameVeryShortVersion() {
        return "LAME3.100";
    }

    getPsyVersion() {
        return "0.93";
    }

    getLameUrl() {
        return "http://www.mp3dev.org/";
    }

    getLameOsBitness() {
        return "32bits";
    }
}

module.exports = Version;
